import sys
import coolfluid as cf
from math import pi

from math import *

u_in = [20.,0.] # 2 dimensioneel
u_wall = [0.,0.]
rho = 1.225
mu = 0.0000178 # p 18

tstep = 5.  #dit is stap in de tijd
num_steps = 50. #aantal stappen

env = cf.Core.environment()
env.log_level = 3

# Basic model setup (container of all the sumilation setup)
model = cf.root.create_component('NavierStokes', 'cf3.solver.ModelUnsteady')

# The domain contains the mesh
domain = model.create_domain()

# Physics contains constants such as viscosity
physics = model.create_physics('cf3.UFEM.KEpsilonPhysics')

# Solver contains all solvers for the problem (conduction, flow, ...)
solver = model.create_solver('cf3.UFEM.Solver')

# Add a concrete Navier-Stokes finite element solver

#add the spalart almares turbulence model
satm = solver.add_unsteady_solver('cf3.UFEM.SpalartAllmaras')

ns_solver = solver.add_unsteady_solver('cf3.UFEM.NavierStokes')

ad_solver = solver.add_unsteady_solver('cf3.UFEM.adjoint.Adjoint')
#ad_solver.ct = Ct0adj
ad_solver.th = 0.4
# ad_solver.area = area
ad_solver.turbulence = 0. # met ke model

gradient1 = solver.add_unsteady_solver('cf3.UFEM.VelocityGradient')
gradient1.options.gradient_tag = 'Adjvelocity_gradient'
gradient1.options.velocity_variable = 'AdjVelocity'
gradient1.options.velocity_tag = 'adjoint_solution'
gradient1.options.gradient_name = 'U'

mesh = domain.load_mesh(file = cf.URI('AlternatieveBuis.msh'), name = 'Mesh')

ns_solver.regions = [mesh.topology.uri()]
ad_solver.regions = [mesh.topology.uri()]
satm.regions = [mesh.topology.uri()]

# initial conditions
solver.InitialConditions.navier_stokes_solution.Velocity = u_in
solver.InitialConditions.density_ratio.density_ratio = 1. # This enables the body force
solver.InitialConditions.adjoint_solution.AdjVelocity = u_in
# viscositeit voor Spalart Allmaras
NU_in = 0.001
solver.InitialConditions.spalart_allmaras_solution.SAViscosity = NU_in


ic_k = solver.InitialConditions.create_initial_condition(builder_name = 'cf3.UFEM.InitialConditionFunction', field_tag = 'ke_solution')#
ic_k.variable_name = 'k'
ic_k.value = [str(1.)]
ic_k.regions = [mesh.topology.uri()]

ic_epsilon = solver.InitialConditions.create_initial_condition(builder_name = 'cf3.UFEM.InitialConditionFunction', field_tag = 'ke_solution')
ic_epsilon.variable_name = 'epsilon'
ic_epsilon.value = [str(1.)]
ic_epsilon.regions = [mesh.topology.uri()]

# adjoint oplossing
ic_k = solver.InitialConditions.create_initial_condition(builder_name = 'cf3.UFEM.InitialConditionFunction', field_tag = 'keAdjoint_solution')
ic_k.variable_name = 'ka'
ic_k.value = [str(1.)]
ic_k.regions = [mesh.topology.uri()]

ic_epsilon = solver.InitialConditions.create_initial_condition(builder_name = 'cf3.UFEM.InitialConditionFunction', field_tag = 'keAdjoint_solution')
ic_epsilon.variable_name = 'epsilona'
ic_epsilon.value = [str(1.)]
ic_epsilon.regions = [mesh.topology.uri()]

# set physical constants
physics.density = rho
physics.dynamic_viscosity = mu

# Boundary conditions
bc = ns_solver.BoundaryConditions
bc.add_constant_bc(region_name = 'inlet', variable_name = 'Velocity').value = u_in
bc.add_constant_bc(region_name = 'outlet', variable_name = 'Pressure').value = 0.
bc.add_constant_bc(region_name = 'top', variable_name = 'Velocity').value = u_wall
bc.add_constant_bc(region_name = 'bottom', variable_name = 'Velocity').value = u_wall
bca = ad_solver.BoundaryConditions

bc_adj_p0 = bca.create_bc_action(region_name = 'outlet', builder_name = 'cf3.UFEM.adjoint.BCAdjointpressurex')
bc_adj_p0.turbulence = 1
bc_adj_p1 = bca.create_bc_action(region_name = 'inlet', builder_name = 'cf3.UFEM.adjoint.BCAdjointpressurex')
bc_adj_p1.turbulence = 1

bc_adj_p3 = bca.create_bc_action(region_name = 'top', builder_name = 'cf3.UFEM.adjoint.BCAdjointpressurey')
bc_adj_p3.turbulence = 1
bc_adj_p4 = bca.create_bc_action(region_name = 'bottom', builder_name = 'cf3.UFEM.adjoint.BCAdjointpressurey')
bc_adj_p4.turbulence = 1

bca.add_constant_bc(region_name = 'outlet', variable_name = 'AdjPressure').value = 0.
bca.add_constant_bc(region_name = 'top', variable_name = 'AdjVelocity').value = u_wall
bca.add_constant_bc(region_name = 'bottom', variable_name = 'AdjVelocity').value = u_wall
bca.add_constant_bc(region_name = 'inlet', variable_name = 'AdjVelocity').value = u_in

NU_wall = 0.

# Boundary conditions for Spalart-Allmaras
	#bc = satm.get_child('BoundaryConditions')
    #bc = satm.BoundaryConditions
	#bc.add_constant_bc(region_name = 'inlet', variable_name = 'SAViscosity').options().set('value', NU_in)
    #bc.add_constant_bc(region_name = 'outlet', variable_name = 'SAViscosity').options().set('value', NU_in)

    #bc.add_constant_bc(region_name = 'bottom', variable_name = 'SAViscosity').options().set('value', NU_wall)
	#bc.add_constant_bc(region_name = 'top', variable_name = 'SAViscosity').options().set('value', NU_wall)

# Solver setup

lss = ns_solver.create_lss(matrix_builder = 'cf3.math.LSS.TrilinosFEVbrMatrix', solution_strategy = 'cf3.math.LSS.TrilinosStratimikosStrategy')
lss.SolutionStrategy.Parameters.preconditioner_type = 'ML'
lss.SolutionStrategy.Parameters.PreconditionerTypes.ML.MLSettings.add_parameter(name = 'ML output', value = 0)
lss.SolutionStrategy.Parameters.PreconditionerTypes.ML.MLSettings.default_values = 'NSSA'
lss.SolutionStrategy.Parameters.PreconditionerTypes.ML.MLSettings.aggregation_type = 'Uncoupled'
lss.SolutionStrategy.Parameters.PreconditionerTypes.ML.MLSettings.smoother_type = 'symmetric block Gauss-Seidel'#'Chebyshev'
lss.SolutionStrategy.Parameters.PreconditionerTypes.ML.MLSettings.smoother_sweeps = 4
lss.SolutionStrategy.Parameters.PreconditionerTypes.ML.MLSettings.smoother_pre_or_post = 'post'

lss.SolutionStrategy.Parameters.LinearSolverTypes.Belos.solver_type = 'Block GMRES'
lss.SolutionStrategy.Parameters.LinearSolverTypes.Belos.SolverTypes.BlockGMRES.convergence_tolerance = 1e-5
lss.SolutionStrategy.Parameters.LinearSolverTypes.Belos.SolverTypes.BlockGMRES.maximum_iterations = 2000
lss.SolutionStrategy.Parameters.LinearSolverTypes.Belos.SolverTypes.BlockGMRES.num_blocks = 1000

lss = ad_solver.LSS
lss.SolutionStrategy.Parameters.preconditioner_type = 'ML'
lss.SolutionStrategy.Parameters.PreconditionerTypes.ML.MLSettings.add_parameter(name = 'ML output', value = 0)
lss.SolutionStrategy.Parameters.PreconditionerTypes.ML.MLSettings.default_values = 'NSSA'
lss.SolutionStrategy.Parameters.PreconditionerTypes.ML.MLSettings.aggregation_type = 'Uncoupled'
lss.SolutionStrategy.Parameters.PreconditionerTypes.ML.MLSettings.smoother_type = 'symmetric block Gauss-Seidel'#'Chebyshev'


lss.SolutionStrategy.Parameters.PreconditionerTypes.ML.MLSettings.smoother_sweeps = 4
lss.SolutionStrategy.Parameters.PreconditionerTypes.ML.MLSettings.smoother_pre_or_post = 'post'

lss.SolutionStrategy.Parameters.LinearSolverTypes.Belos.solver_type = 'Block GMRES'
lss.SolutionStrategy.Parameters.LinearSolverTypes.Belos.SolverTypes.BlockGMRES.convergence_tolerance = 1e-5
lss.SolutionStrategy.Parameters.LinearSolverTypes.Belos.SolverTypes.BlockGMRES.maximum_iterations = 2000
lss.SolutionStrategy.Parameters.LinearSolverTypes.Belos.SolverTypes.BlockGMRES.num_blocks = 1000


# time setup
time = model.create_time()
time.time_step = tstep
time.end_time = num_steps*tstep

# run the simulation (forward, NS only)
#solver.TimeLoop.options.disabled_actions = ['Adjoint']
model.simulate()

domain.write_mesh(cf.URI('out-before-adjoint.pvtu'))

# run adjoint, starting from converged NS solution
#solver.TimeLoop.options.disabled_actions = ['NavierStokes'] # NS disabled, adjoint enabled
#solver.options.disabled_actions = ['InitialConditions'] # disable initial conditions
#time.time_step = tstep/10.
#time.end_time += num_steps*tstep #eindtijd
#model.simulate()

#domain.write_mesh(cf.URI('out-after-adjoint3.pvtu'))

# writer = domain.create_component('VTKWriter', 'cf3.vtk.MultiblockWriter')
# writer.mesh = mesh
# writer.file = cf.URI('outputadjoint2initlam2-{e}.vtm'.format(e=u_in[0]))
# writer.execute()
