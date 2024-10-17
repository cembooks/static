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

#include "deal.II/base/multithread_info.h"

#include <iostream>
#include <string>

#include "misc.hpp"
#include "solver.hpp"

using namespace Misc;

/**
 * \brief This is a wrap-around class. It contains the main loop of the program
 * that implements the
 * [Floating conductor (flc/)](@ref page_flc)
 * numerical experiment.
 *****************************************************************************/
class BatchFLC : public SettingsFLC
{
public:
	void run()
	{
		if (nr_threads_max > 0)
			MultithreadInfo::set_thread_limit(nr_threads_max);

		Assert( DIMENSION__ > 1 , ExcInternalError());
		Assert( DIMENSION__ < 4 , ExcInternalError());

		std::string dir;

		dir = (DIMENSION__ == 2) ? "Data/ring/" : "Data/shell/";

		std::cout
			<< "Program: flc\n"
			<< "Dimensions: " << DIMENSION__ << "\n"
			<< "Writing to: " << dir << "\n";

		MainOutputTable table(DIMENSION__);

		for (unsigned int p = 1; p < 4; p++)
		{
			table.clear();
#if DIMENSION__ == 2
			for (unsigned int r = 2; r < 6; r++)
#endif
#if DIMENSION__ == 3
			for (unsigned int r = 3; r < 7; r++)
#endif
			{
				table.add_value("r", r);
				table.add_value("p", p);

				SolverFLC<DIMENSION__> problem(p, 2, r,
					dir + "solution_p" + std::to_string(p) +
					"_r" + std::to_string(r));
				table.add_value("ndofs", problem.get_n_dofs());
				table.add_value("ncells", problem.get_n_cells());
				table.add_value("L2", problem.get_L2_norm());
				table.add_value("H1", problem.get_H1_norm());
			}
			table.save(dir + "table_PHI_p" + std::to_string(p));
		}
	}
};

int main()
{
	try
	{
		BatchFLC batch;
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

