/******************************************************************************/
/*           GOLEM - Multiphysics of faulted geothermal reservoirs            */
/*                                                                            */
/*          Copyright (C) 2017 by Antoine B. Jacquey and Mauro Cacace         */
/*             GFZ Potsdam, German Research Centre for Geosciences            */
/*                                                                            */
/*    This program is free software: you can redistribute it and/or modify    */
/*    it under the terms of the GNU General Public License as published by    */
/*      the Free Software Foundation, either version 3 of the License, or     */
/*                     (at your option) any later version.                    */
/*                                                                            */
/*       This program is distributed in the hope that it will be useful,      */
/*       but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/*        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       */
/*                GNU General Public License for more details.                */
/*                                                                            */
/*      You should have received a copy of the GNU General Public License     */
/*    along with this program.  If not, see <http://www.gnu.org/licenses/>    */
/******************************************************************************/

#include "GolemConditionalFunctionEnableControl.h"
#include "Function.h"

registerMooseObject("GolemApp", ConditionalFunctionEnableControl);

InputParameters
ConditionalFunctionEnableControl::validParams()
{
  InputParameters params = ConditionalEnableControl::validParams();

  params.addRequiredParam<FunctionName>("conditional_function",
                                        "The function to give a true or false value");

  params.addClassDescription(
      "Control for enabling/disabling objects when a function value is true");

  return params;
}

ConditionalFunctionEnableControl::ConditionalFunctionEnableControl(
    const InputParameters & parameters)
  : ConditionalEnableControl(parameters), _function(getFunction("conditional_function"))
{
}

bool
ConditionalFunctionEnableControl::conditionMet(const unsigned int & /*i*/)
{
  return _function.value(_t);
}
