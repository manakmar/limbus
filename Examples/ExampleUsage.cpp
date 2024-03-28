// LIMBUS 1.0 - the Library for Measuring Ball Union Structures
// ------------------------------------------------------------
// This program computes the volume of the union of a set of balls.
// The volume is computed analytically, using the regular triangulation
// and alpha-shapes functionality from the CGAL library. Analytic equations
// come from the paper "Measuring Space Filling Diagrams and Voids" 
// written by Herbert Edelsbrunner and Ping Fu in 1994 (available at 
// https://pub.ista.ac.at/~edels/Papers/1994-07-MeasuringSpaceFillingDiagrams.pdf).
// Alternatively, the volume can also be roughly estimated by Monte Carlo 
// sampling.
//
// Prerequisites (Fedora 39 Linux)
// -------------------------------
// CGAL-devel-5.6.1-1
// boost-1.81.0-8
// boost-program-options-1.81.0-8
// boost-devel-1.81.0-8
//
// sudo dnf install CGAL-devel boost boost-program-options boost-devel
//
// Compilation
// -----------
// g++ -I ../ ExampleUsage.cpp -o limbus -lgmp -lgmpxx -lboost_program_options
//
// Usage
// -----
// ./limbus -a -m < input_xyzr.txt  

#include <cstdlib>
#include <iostream>
#include <vector>

#include <boost/program_options.hpp>

#include <limbus/Limbus.h>

using namespace limbus;
using namespace std;

void readBalls(istream& in, vector<Ball>& out_balls);

int main(int argc, char *argv[])
{
  namespace po = boost::program_options;
  
  ios::sync_with_stdio(false); // faster
  
  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
      ("help,h", "print a help message")
      ("analytic,a", "analytically computed volume via regular triangulation and alpha shapes")
      ("sampling,s", "numerically estimated volume via Monte Carlo sampling")
  ;
  
  // Parse program options
  po::variables_map vm;
  po::store(
    po::command_line_parser(argc,argv)
        .options(desc)
        .style(
            po::command_line_style::unix_style | 
            po::command_line_style::allow_short |
            po::command_line_style::allow_long |
            po::command_line_style::allow_dash_for_short)
        .run(),
    vm);
  vm.notify();  

  if (vm.count("help") >= 1)
  {
    cerr << desc << endl;
    return EXIT_SUCCESS;
  }
    
  vector<Ball> balls;
  readBalls(cin, balls);
  
  if (vm.count("analytic") >= 1 || vm.count("sampling") == 0)
  {
    UnionBalls ub(balls.begin(), balls.end());
    cout << "Analytically computed volume: " << ub.computeVolume() << endl;
  }
  if (vm.count("sampling") >= 1)
  {
    UnionBallsMonteCarlo mc(balls.data(), balls.size());
    cout << "Monte Carlo sampling volume:  " << mc.computeVolume() << endl;
  }
  
  return EXIT_SUCCESS;
}

void readBalls(istream& in, vector<Ball>& out_balls)
{
  string line;
  while (getline(cin, line))
  {
    // Skip lines without content
    bool is_white_line = true;
    for (char c: line)
    {
      is_white_line &= std::isspace(c);
    }
    if (is_white_line) {
      continue;
    }
    
    double x, y, z, r;
    std::istringstream(line) >> std::skipws >> x >> y >> z >> r;
    out_balls.push_back(Ball(x, y, z, r));
  }
}
