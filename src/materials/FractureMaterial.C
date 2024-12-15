#include "FractureMaterial.h"
#include "Assembly.h"
#include "MooseMesh.h"
#include "libmesh/quadrature.h"
#include "Function.h"
#include "GolemH.h"

registerMooseObject("GolemApp", FractureMaterial);

InputParameters
FractureMaterial::validParams()
{
  InputParameters params = GolemMaterialBase::validParams();
  params.addRequiredCoupledVar("displacements", "The displacement variables vector.");
  params.addCoupledVar("pore_pressure", "The pore pressure variable.");
  params.addParam<RealVectorValue>("body_force", RealVectorValue(0, 0, 0), "The body force vector F_b [N/m^3].");
  params.addParam<bool>("include_surface_terms", true, "Flag to include boundary contributions.");
  params.addParam<Real>("normal_stress_closure", 1e6, "Normal effective stress \sigma'_o causing fracture closure.");
  params.addParam<Real>("zero_contact_aperture", 1e-3, "Aperture af_0 at zero contact stress.");
  params.addParam<Real>("youngs_modulus", 1e10, "Young's modulus E of the rock.");
  params.addParam<Real>("mesh_element_size", 1.0, "Mesh element size h.");
  return params;
}

FractureMaterial::FractureMaterial(const InputParameters & parameters)
  : GolemMaterialBase(parameters),
    _ndisp(coupledComponents("displacements")),
    _disp(3),
    _grad_disp(3),
    _stress(declareProperty<RankTwoTensor>("stress")),
    _body_force(getParam<RealVectorValue>("body_force")),
    _include_surface_terms(getParam<bool>("include_surface_terms")),
    _sigma_o(getParam<Real>("normal_stress_closure")),
    _af_0(getParam<Real>("zero_contact_aperture")),
    _youngs_modulus(getParam<Real>("youngs_modulus")),
    _mesh_element_size(getParam<Real>("mesh_element_size"))
{
  if (_ndisp != _mesh.dimension())
    mooseError("The number of displacement variables must match the mesh dimension!");
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    _disp[i] = &coupledValue("displacements", i);
    _grad_disp[i] = &coupledGradient("displacements", i);
  }
}

void FractureMaterial::computeQpProperties()
{
  computeStrain();
  computeStress();
  assembleResidual();
}

void FractureMaterial::computeStrain()
{
  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    RankTwoTensor grad_u = RankTwoTensor::initializeFromRows(
      (*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp], (*_grad_disp[2])[_qp]);
    _strain[_qp] = 0.5 * (grad_u + grad_u.transpose());
  }
}

void FractureMaterial::computeStress()
{
  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    _stress[_qp] = _Cijkl[_qp] * _strain[_qp];
  }
}

void FractureMaterial::assembleResidual()
{
  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    RealVectorValue residual = RealVectorValue(0, 0, 0);

    // Body force term
    residual -= _body_force * _q_point[_qp];

    // Internal stress contribution
    for (unsigned int i = 0; i < _ndisp; ++i)
      for (unsigned int j = 0; j < _ndisp; ++j)
        residual[i] += _stress[_qp](i, j) * (*_grad_disp[j])[_qp];

    // Fracture mechanics: normal stiffness contribution
    Real jump = computeDisplacementJump();
    if (_af_0 >= jump)
    {
      Real normal_stiffness = (_af_0 * _sigma_o) / (9 * std::pow(jump, 2));
      residual -= normal_stiffness * jump * _q_point[_qp];
    }

    // Boundary contributions
    if (_include_surface_terms)
    {
      for (unsigned int side = 0; side < _mesh.n_boundaries(); ++side)
      {
        const BoundaryInfo & binfo = _mesh.get_boundary_info(side);
        for (auto & face : binfo)
        {
          // Neumann contributions from sigma \hat{n} and pf \hat{n}
          RealVectorValue boundary_stress = _stress[_qp] * face.normal();
          if (isFluidBoundary(face))
            boundary_stress += coupledValue("pore_pressure") * face.normal();

          residual += boundary_stress * _JxW[_qp];
        }
      }
    }

    _residual[_qp] = residual;
  }
}

Real FractureMaterial::computeDisplacementJump()
{
  RealVectorValue jump = (*_disp[0])[_qp] - (*_disp[1])[_qp];
  return jump * _normal[_qp];
}

bool FractureMaterial::isFluidBoundary(const BoundaryInfo & binfo) const
{
  return binfo.boundary_id() == _fluid_boundary_id;
}
