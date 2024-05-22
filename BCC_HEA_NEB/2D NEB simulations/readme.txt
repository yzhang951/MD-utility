input.lmp:			NEB script

build.lmp & build_final.lmp:	Energy minimization for initial and final states: 
(First use 'build.lmp' to relieve stress for both initial and final states,
then use 'build_final.lmp' to keep the final state box the same size as the initial, since 'fix box/relax tri' is used previously.)

calc 1-3: initial and final configuration files after minimization and the resulting .cfg files

