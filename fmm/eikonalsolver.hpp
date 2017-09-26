/*! \class EikonalSolver
    \brief Abstract class that serves as interface for the actual EikonalSolvers implemented.
    It requires (at least) the computeInternal method to be implemented,

    It uses as a main container the nDGridMap class. The nDGridMap template paramenter
    has to be an FMCell or something inherited from it.

    Copyright (C) 2015 Javier V. Gomez
    www.javiervgomez.com

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef EIKONALSOLVER_H_
#define EIKONALSOLVER_H_

#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <array>
#include <chrono>

#include <boost/concept_check.hpp>

#include "solver.hpp"
#include "../console/console.h"


#include <algorithm>
#include <iterator>

template <class grid_t>
class EikonalSolver : public Solver<grid_t>{

    public:
        EikonalSolver() : Solver<grid_t>("EikonalSolver") {}
        EikonalSolver(const std::string& name) : Solver<grid_t>(name) {}

        /** \brief Solves nD Eikonal equation for cell idx. If heuristics are activated, it will add
            the estimated travel time to goal with current velocity. */
        /*const */double/*&*/ solveEikonal
        (const int & idx) {
	    
	    const auto& cellI        = grid_->getCell(idx);
            const auto& cellIArrTime = cellI.getArrivalTime();
            double* TvalIter_ = &Tvalues_[0];
            size_t a = 0;// a parameter of the Eikonal equation.
            for ( size_t dim = 0; dim < grid_t::getNDims(); ++dim, ++TvalIter_ ) {
                const double& minTInDim = grid_->getMinValueInDim(idx, dim);
                *TvalIter_ = minTInDim;
                if ( minTInDim < cellIArrTime ) { a++; }
            }
            
	    
            // Sort the neighbor values to make easy the following code.
            /// \todo given that this sorts a small vector, a n^2 methods could be better. Test it.
            std::sort(&Tvalues_[0], &Tvalues_[grid_t::getNDims()]);
	    
	    double updatedT_ = std::numeric_limits<double>::infinity();
            if ( a > 0 ) {
                scaledT_ = grid_->getLeafSize() / cellI.getVelocity();
                TvalIter_ = &Tvalues_[0];
                sumTT_ = (*TvalIter_)*(*TvalIter_);
                sumT_  = (*TvalIter_++);
                /*= */solveEikonal1Dim(updatedT_/*cellI*/);
//                 if ( updatedT_ - (*TvalIter_) < utils::COMP_MARGIN ) { return updatedT_; } // If increasing one dimension will not improve time.
                scaledTSqr_ = scaledT_*scaledT_;
                for (size_t i = 2; i <= a; ++i, ++TvalIter_) {
                    sumTT_ += (*TvalIter_)*(*TvalIter_);
                    sumT_  += (*TvalIter_);
                    /*updatedT_ = */solveEikonalNDims(updatedT_,/*cellI, */i);
                    if ( updatedT_ - (*TvalIter_) < utils::COMP_MARGIN ) { break; } // If increasing one dimension will not improve time.
                }
            }
            return updatedT_;
        }

        
    protected:
        /** \brief Solves the Eikonal equation assuming that Tvalues_
            is sorted. */
//         template<typename TCellType>
        /*double */void solveEikonal1Dim
        (/*const TCellType& cellI*/double& updT) const {
            // Solve for 1 dimension.
            updT = Tvalues_[0] + /*grid_->getLeafSize() / cellI.getVelocity()*/scaledT_;
        }
// 	template<typename TCellType>
        /*double */void solveEikonalNDims
        (/*const TCellType& cellI,*/double& updT, const size_t& dim) const {
            // These a,b,c values are simplified since leafsize^2, which should be present in the three
            // terms but they are cancelled out when solving the quadratic function.
            const double&& aTime2 = 2*dim;
            const double&& minB   = 2*sumT_;
            const double&& c      = sumTT_ - scaledTSqr_;
            const double&& quad_term = minB*minB - 2*aTime2*c;

            updT = (quad_term >= 0) ? (minB + sqrt(quad_term))/aTime2 : std::numeric_limits<double>::infinity();
            // for infinity grid_->getLeafSize() * grid_->getLeafSize()/(grid_->getCell(idx).getVelocity()*grid_->getCell(idx).getVelocity()) + Tvalues_[0];
        }

        /** \brief Auxiliar vector with values T0,T1...Tn-1 variables in the Discretized Eikonal Equation. */
        std::array<double, grid_t::getNDims()> Tvalues_;
//         double* TvalIter_;
	double sumT_;
	double sumTT_;
// 	double updatedT_;
        double scaledT_;
        double scaledTSqr_;

        /** \brief Auxiliar array which stores the neighbor of each iteration of the computeFM() function. */
//         std::array <unsigned int, 2*grid_t::getNDims()> neighbors_;

        using Solver<grid_t>::grid_;
};

#endif /* EIKONALSOLVER_H_*/
