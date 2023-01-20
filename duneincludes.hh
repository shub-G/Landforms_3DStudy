#include<dune/common/parametertreeparser.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parallel/mpiguard.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/timer.hh>

#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/vtk/vtksequencewriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/utility/structuredgridfactory.hh>

#if HAVE_ALBERTA
#include<dune/grid/albertagrid.hh>
#include<dune/grid/albertagrid/dgfparser.hh>
#endif
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#if HAVE_ALUGRID
#include<dune/grid/alugrid.hh>
#include<dune/grid/io/file/dgfparser/dgfalu.hh>
#include<dune/grid/io/file/dgfparser/dgfparser.hh>
#endif

#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/superlu.hh>

#include<dune/pdelab.hh>

//#include<dune/pdelab/newton/newton.hh>
//#include<dune/pdelab/common/function.hh>
//#include<dune/pdelab/common/vtkexport.hh>
//#include<dune/pdelab/common/geometrywrapper.hh>
//#include<dune/pdelab/common/quadraturerules.hh>
//
//#include<dune/pdelab/finiteelement/localbasiscache.hh>
//#include<dune/pdelab/finiteelementmap/p0fem.hh>
//#include<dune/pdelab/finiteelementmap/qkdg.hh>
//
//#include<dune/pdelab/constraints/conforming.hh>
//#include<dune/pdelab/constraints/noconstraints.hh>
//#include<dune/pdelab/constraints/p0.hh>
//#include<dune/pdelab/constraints/common/constraints.hh>
//#include<dune/pdelab/constraints/common/constraintsparameters.hh>
//#include<dune/pdelab/finiteelementmap/qkfem.hh>
//#include<dune/pdelab/function/callableadapter.hh>
////
//#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
//#include<dune/pdelab/gridfunctionspace/subspace.hh>
//#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
//#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
//#include<dune/pdelab/gridfunctionspace/interpolate.hh>
//#include<dune/pdelab/gridfunctionspace/localvector.hh>
//#include<dune/pdelab/gridfunctionspace/vtk.hh>
//
//#include<dune/pdelab/gridoperator/gridoperator.hh>
//#include<dune/pdelab/gridoperator/onestep.hh>
//
//#include<dune/pdelab/backend/istl.hh>
//
//#include<dune/pdelab/stationary/linearproblem.hh>
//#include<dune/pdelab/instationary/onestep.hh>
//
//#include<dune/pdelab/localoperator/defaultimp.hh>
//#include<dune/pdelab/localoperator/pattern.hh>
//#include<dune/pdelab/localoperator/flags.hh>
//#include<dune/pdelab/localoperator/idefault.hh>

#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
