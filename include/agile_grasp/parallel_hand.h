/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2015, Andreas ten Pas
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef PARALLEL_HAND_H
#define PARALLEL_HAND_H

#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "finger_hand.h"

/** ParallelHand class
 *
 * \brief Calculate collision-free fingers
 * 
 * This class calculates a set of collision-free finger placements. The parameters are the dimensions of the 
 * robot hand and the desired "bite" that the grasp must have. The bite is how far the hand can be 
 * moved into the object.
 * 
*/
class ParallelHand : public FingerHand
{
public:
	/**
	* \brief Default constructor.
	*/
	ParallelHand()
	{
	}

	/**
	 * \brief Constructor.
	 * \param finger_width the width of the fingers
	 * \param hand_outer_diameter the maximum aperture of the robot hand
	 * \param hand_depth the length of the fingers
	*/
	ParallelHand(double finger_width, double hand_outer_diameter, double hand_depth);

	/**
	* \brief Find collision-free finger placements.
	* \param bite the minimum object height
	*/
	void evaluateFingers(double bite);

	/**
	* \brief Find robot hand configurations that fit the cloud.
	*/
	void evaluateHand();

	/**
	* \brief Set the grasp parameters.
	*
	* The parameter @p bite is used to calculate the grasp width by only evaluating the width of 
	* the points below @p bite.
	*
	* \param bite the minimum object height
	*/
	void evaluateGraspParameters(double bite);

	/**
	* \brief Select center hand and try to extend it into the object as deeply as possible.
	* \param init_deepness the initial deepness (usually the minimum object height)
	* \param max_deepness the maximum allowed deepness (usually the finger length)
	*/
	void deepenHand(double init_deepness, double max_deepness);

};

#endif
