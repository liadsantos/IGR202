// ----------------------------------------------------------------------------
// RigidSolver.hpp
//
//  Created on: 18 Dec 2020
//      Author: Kiwon Um
//        Mail: kiwon.um@telecom-paris.fr
//
// Description: Simple Rigid Body Solver (DO NOT DISTRIBUTE!)
//
// Copyright 2020-2023 Kiwon Um
//
// The copyright to the computer program(s) herein is the property of Kiwon Um,
// Telecom Paris, France. The program(s) may be used and/or copied only with
// the written permission of Kiwon Um or in accordance with the terms and
// conditions stipulated in the agreement/contract under which the program(s)
// have been supplied.
// ----------------------------------------------------------------------------

#ifndef _RIGIDSOLVER_HPP_
#define _RIGIDSOLVER_HPP_

#include <glm/ext/matrix_transform.hpp>

#include "Vector3.hpp"
#include "Matrix3x3.hpp"
#include "quaternion.hpp"

struct BodyAttributes {
  BodyAttributes() :
    X(0, 0, 0), R(Mat3f::I()), P(0, 0, 0), L(0, 0, 0), Q(Quaternion(1, 0, 0, 0)),
    V(0, 0, 0), omega(0, 0, 0), omegaTilde(Quaternion(1, 0, 0, 0)),
    F(0, 0, 0), tau(0, 0, 0) {}

  glm::mat4 worldMat() const
  {
    return glm::mat4(           // column-major
      R(0,0), R(1,0), R(2,0), 0,
      R(0,1), R(1,1), R(2,1), 0,
      R(0,2), R(1,2), R(2,2), 0,
      X[0],   X[1],   X[2],   1);
  }

  tReal M;                      // mass
  Mat3f I0, I0inv;              // inertia tensor and its inverse in body space
  Mat3f Iinv;                   // inverse of inertia tensor

  // rigid body state
  Vec3f X;                      // position
  Mat3f R;                      // rotation
  Vec3f P;                      // linear momentum
  Vec3f L;                      // angular momentum
  Quaternion Q;                 // quaternion

  // auxiliary quantities
  Vec3f V;                      // linear velocity
  Vec3f omega;                  // angular velocity
  Quaternion omegaTilde;        // angular velocity for quaternion

  // force and torque
  Vec3f F;                      // force
  Vec3f tau;                    // torque

  // mesh's vertices in body space
  std::vector<Vec3f> vdata0;
};

class Box : public BodyAttributes {
public:
  explicit Box(
    const tReal w=1.0, const tReal h=1.0, const tReal d=1.0, const tReal dens=10.0,
    const Vec3f v0=Vec3f(0, 0, 0), const Vec3f omega0=Vec3f(0, 0, 0)) :
    width(w), height(h), depth(d)
  {
    V = v0;                     // initial velocity
    omega = omega0;             // initial angular velocity

    // Calculate physical attributes
    M = dens * w * h * d;       // mass = density * volume
    I0 = Mat3f(Vec3f(           // inertia tensor on body space
      1.f/12.f * M * (pow(h,2) + pow(d,2)),
      1.f/12.f * M * (pow(w,2) + pow(d,2)),
      1.f/12.f * M * (pow(w,2) + pow(h,2)))
    );
    I0inv = I0.inverse();       // inverse of inertia tensor on body space
    Iinv = R * I0inv * R.transpose();

    // vertices data (8 vertices)
    vdata0.push_back(Vec3f(-0.5*w, -0.5*h, -0.5*d));
    vdata0.push_back(Vec3f( 0.5*w, -0.5*h, -0.5*d));
    vdata0.push_back(Vec3f( 0.5*w,  0.5*h, -0.5*d));
    vdata0.push_back(Vec3f(-0.5*w,  0.5*h, -0.5*d));

    vdata0.push_back(Vec3f(-0.5*w, -0.5*h,  0.5*d));
    vdata0.push_back(Vec3f( 0.5*w, -0.5*h,  0.5*d));
    vdata0.push_back(Vec3f( 0.5*w,  0.5*h,  0.5*d));
    vdata0.push_back(Vec3f(-0.5*w,  0.5*h,  0.5*d));
  }

  // rigid body property
  tReal width, height, depth;
};

class RigidSolver {
public:
  explicit RigidSolver(
    BodyAttributes *body0=nullptr, const Vec3f g=Vec3f(0, 0, 0)) :
    body(body0), _g(g), _step(0), _sim_t(0) {}

  void init(BodyAttributes *body0)
  {
    body = body0;
    _step = 0;
    _sim_t = 0;
  }

  void step(const tReal dt)
  {
    std::cout << "t=" << _sim_t << " (dt=" << dt << ")" << std::endl;

    computeForceAndTorque();

    /// Time integration - temporal evolution of state variables
    // linear
    body->X = body->X + dt * body->P / body->M;     // position
    body->P = body->P + dt * body->F;               // linear momentum
    body->L = body->L + dt * body->tau;             // angular momentum

    // angular
    body->Iinv = body->R * body->I0inv * body->R.transpose();
    body->omega = body->Iinv * body->L;             // angular velocity
    //body->R = body->R + dt * (body->omega).crossProductMatrix() * body->R;  // rotation matrix

    // quaternion
    body->omegaTilde = Quaternion(body->omega, 0);
    body->Q = body->Q + Quaternion(body->omega, 0).product(body->Q) * dt / 2.f;
    body->R = body->Q.normalize().rotationMatrix();

    ++_step;
    _sim_t += dt;
  }

  BodyAttributes *body;

private:
  void computeForceAndTorque()
  {
    // Force (mass * accelaration)
    body->F = body->M * _g;
    body->tau = Vec3f(0, 0, 0);

    // Instance force at the very first step
    if(_step == 1) {
      body->F += Vec3f(2.15, 0.25, 0.03);

      // torque (r - x) crossProduct F
      body->tau = (body->R * body->vdata0[0]).crossProduct(body->F); // apply only when the force is applied
    }
  }

  // simulation parameters
  Vec3f _g;                     // gravity
  tIndex _step;                 // step count
  tReal _sim_t;                 // simulation time
};

#endif  /* _RIGIDSOLVER_HPP_ */
