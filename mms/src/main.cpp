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
#include "project_PHI_to_E.hpp"
#include "project_PHI_to_D.hpp"

using namespace Misc;

/**
 * \brief This is a wrap-around class. It contains the main loop of the program
 * that implements the
 * [Method of manufactured solutions (mms/)](@ref page_mms)
 * numerical experiment.
 *****************************************************************************/
class BatchMMS : public SettingsMMS
{
public:
	void run()
	{
		if (nr_threads_max > 0)
			MultithreadInfo::set_thread_limit(nr_threads_max);

		Assert( DIMENSION__ > 1 , ExcInternalError())
		Assert( DIMENSION__ < 4 , ExcInternalError())

		std::string dir;
		std::string fname;

		if (HYPERCUBE__ == 1)
		{
			dir = ( DIMENSION__ == 2 ) ? "Data/square/" : "Data/cube/";
		}else
		{
			dir = ( DIMENSION__ == 2 ) ? "Data/circle/" : "Data/sphere/";
		}

		std::cout
			<< "Program: mms\n"
			<< "Dimensions: " << DIMENSION__ << "\n"
			<< "Writing to: " << dir << "\n";

		MainOutputTable table_PHI(DIMENSION__);
		MainOutputTable table_E(DIMENSION__);
		MainOutputTable table_D(DIMENSION__);

		for (unsigned int p = 1; p < 4; p++)
		{
			table_PHI.clear();
			table_E.clear();
			table_D.clear();

#if DIMENSION__ == 2
		for (unsigned int r = 9; r < 13; r++)
#endif
#if DIMENSION__ == 3
		for (unsigned int r = 8; r < 12; r++)
#endif
			{
				fname = dir + "solution_PHI_p" + std::to_string(p)
					+ "_r" + std::to_string(r);

				// Calculating potential
				table_PHI.add_value("r", r);
				table_PHI.add_value("p", p);

				if (SettingsMMS::print_time_tables)
					std::cout << "Time table PHI \n";

				SolverMMS<DIMENSION__> problem(p, 1, r, fname);
				table_PHI.add_value("ndofs", problem.get_n_dofs());
				table_PHI.add_value("ncells", problem.get_n_cells());
				table_PHI.add_value("L2", problem.get_L2_norm());
				table_PHI.add_value("H1", problem.get_H1_norm());

				problem.clear();

				{ // Calculating electrostatic field
					fname = dir + "solution_E_p" + std::to_string(p)
						+ "_r" + std::to_string(r);

					table_E.add_value("r", r);
					table_E.add_value("p", p);

					if (SettingsMMS::print_time_tables)
						std::cout << "Time table E \n";

					ExactSolutionMMS_E<DIMENSION__> exact_solution;

					ProjectPHItoE projector(
						p-1,
						problem.get_mapping_degree(),
						problem.get_tria(),
						problem.get_dof_handler(),
						problem.get_solution(),
						fname,
						& exact_solution,
						false,
						Settings::print_time_tables,
						Settings::project_exact_solution,
						Settings::log_cg_convergence);

					table_E.add_value("ndofs", projector.get_n_dofs());
					table_E.add_value("ncells", projector.get_n_cells());
					table_E.add_value("L2", projector.get_L2_norm());
					table_E.add_value("H1", 0.0);
				}
				{ // Calculating displacement
					fname = dir + "solution_D_p" + std::to_string(p)
						+ "_r" + std::to_string(r);

					table_D.add_value("r", r);
					table_D.add_value("p", p);

					if (SettingsMMS::print_time_tables)
						std::cout << "Time table D \n";

					ExactSolutionMMS_D<DIMENSION__> exact_solution;

					ProjectPHItoD projector(
						p-1,
						problem.get_mapping_degree(),
						problem.get_tria(),
						problem.get_dof_handler(),
						problem.get_solution(),
						fname,
						& exact_solution,
						false,
						Settings::print_time_tables,
						Settings::project_exact_solution,
						Settings::log_cg_convergence);

					table_D.add_value("ndofs", projector.get_n_dofs());
					table_D.add_value("ncells", projector.get_n_cells());
					table_D.add_value("L2", projector.get_L2_norm() / ep_0);
					table_D.add_value("H1", 0.0);
				}
			}
			// Saving convergence tables
			std::cout << "Table PHI\n";
			table_PHI.save(dir + "table_PHI_p" + std::to_string(p));

			std::cout << "Table E\n";
			table_E.save(dir + "table_E_p" + std::to_string(p));

			std::cout << "Table D\n";
			table_D.save(dir + "table_D_p" + std::to_string(p));
		}
	}
};

int main()
{
	try
	{
		BatchMMS batch;
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

