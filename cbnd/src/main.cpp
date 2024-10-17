/******************************************************************************
* Copyright (C) Siarhei Uzunbajakau, 2023.
*
* This program is free software. You can use, modify, and redistribute it under
* the terms of the GNU Lesser General Public License as published by the Free
* Software Foundation, either version 3 or (at your option) any later version.
* This program is distributed without any warranty.
*
* Refer to COPYING.LESSER for more details.
******************************************************************************/

#define BOOST_ALLOW_DEPRECATED_HEADERS

#include <deal.II/base/multithread_info.h>

#include <iostream>
#include <string>

#include "misc.hpp"
#include "solver.hpp"

using namespace Misc;

/**
 * \brief This is a wrap-around class. It contains the main loop of the program
 * that implements the
 * [Effect of curved boundaries (cbnd/)](@ref page_cbnd)
 * numerical experiment.
 *****************************************************************************/
class BatchCBND : public SettingsCBND
{
public:
	void run()
	{
		if (nr_threads_max > 0)
			MultithreadInfo::set_thread_limit(nr_threads_max);

		Assert( DIMENSION__ > 1 , ExcInternalError());
		Assert( DIMENSION__ < 4 , ExcInternalError());

		std::string dir;

#if DIMENSION__ == 2
		dir = (IS_BC_EXACT__ == 1) ? "Data/ring-exact/" : "Data/ring/";
#endif

#if DIMENSION__ == 3
		dir = (IS_BC_EXACT__ == 1) ? "Data/shell-exact/" : "Data/shell/";
#endif

		std::cout
			<< "Program: cbnd\n"
			<< "Dimensions: " << DIMENSION__ << "\n"
			<< "Condition: " << IS_BC_EXACT__ << "\n"
			<< "Writing to: " << dir << "\n";

		MainOutputTable table_PHI(DIMENSION__);

		for (unsigned int p = 1; p < 4; p++)
		{
			table_PHI.clear();

#if DIMENSION__ == 2
			for (unsigned int r = 15; r < 19; r++)
#endif
#if DIMENSION__ == 3
			for (unsigned int r = 9; r < 13; r++)
#endif
			{ // Calculating potential
				table_PHI.add_value("r", r);
				table_PHI.add_value("p", p);

//				SolverCBND<DIMENSION__> problem(p, 2, r,
				SolverCBND<DIMENSION__> problem(p, 1, r,
					dir + "solution_PHI_p" + std::to_string(p) +
					"_r" + std::to_string(r));
				table_PHI.add_value("ndofs", problem.get_n_dofs());
				table_PHI.add_value("ncells", problem.n_cells);
				table_PHI.add_value("L2", problem.get_L2_norm());
				table_PHI.add_value("H1", problem.get_H1_norm());

				problem.clear();
			}
			// Saving convergence tables
			std::cout << "Table PHI\n";
			table_PHI.save(dir + "table_PHI_p" + std::to_string(p));
		}
	}
};

int main()
{
	try
	{
		BatchCBND batch;
		batch.run();
	}
	catch (std::exception &exc)
	{
		std::cerr
			<< std::endl
			<< std::endl
			<< "----------------------------------------------------"
			<< std::endl;
		std::cerr
			<< "Exception on processing: " << std::endl
			<< exc.what() << std::endl
			<< "Aborting!" << std::endl
			<< "----------------------------------------------------"
			<< std::endl;
		return 1;
	}
	catch (...)
	{
		std::cerr
			<< std::endl
			<< std::endl
			<< "----------------------------------------------------"
			<< std::endl;
		std::cerr
			<< "Unknown exception!" << std::endl
			<< "Aborting!" << std::endl
			<< "----------------------------------------------------"
			<< std::endl;
		return 1;
	}

	return 0;
}

