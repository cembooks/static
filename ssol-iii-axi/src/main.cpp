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
#include "project_Aphi_to_Hrz.hpp"
#include "project_Aphi_to_Brz.hpp"

using namespace Misc;

/**
 * \brief This is a wrap-around class. It contains the main loop of the program
 * that implements the
 * [Axisymmetric - thick spherical coil with magnetic core (ssol-iii-axi/)](@ref page_ssol_iii_axi)
 * numerical experiment.
 *****************************************************************************/
class BatchSSOLIIIAXI : public SettingsSSOLIIIAXI
{
public:
	void run()
	{
		if (nr_threads_max > 0)
			MultithreadInfo::set_thread_limit(nr_threads_max);

		std::string dir = "Data/circle/";
		std::string fname;

		std::cout
			<< "Program: ssol-iii-axi\n"
			<< "Dimensions: " << 2 << "\n"
			<< "Writing to: " << dir << "\n";

		MainOutputTable table_H(2);
		MainOutputTable table_B(2);

		for (unsigned int p = 1; p < 4; p++)
		{
			table_H.clear();
			table_B.clear();

			for (unsigned int r = 15; r < 19; r++)
			{ // Calculate A'.
				fname = dir + "solution_A_p" + std::to_string(p) +
					"_r" + std::to_string(r);


				if (SettingsSSOLIIIAXI::print_time_tables)
					std::cout << "Time table A \n";

				SolverSSOLIIIAXI problem(p, 2, r, fname);

				problem.clear();
				{ // Calculate H'.
					fname = dir + "solution_H_p" + std::to_string(p) +
						"_r" + std::to_string(r);

					table_H.add_value("r", r);
					table_H.add_value("p", p);

					if (SettingsSSOLIIIAXI::print_time_tables)
						std::cout << "Time table H \n";

					ExactSolutionSSOLIIIAXI_H exact_solution;

					ProjectAphiToHrz projector(
						p-1,
						problem.get_mapping_degree(),
						problem.get_tria(),
						problem.get_dof_handler(),
						problem.get_solution(),
						fname,
						& exact_solution,
						Settings::print_time_tables,
						Settings::project_exact_solution,
						Settings::log_cg_convergence);

					table_H.add_value("ndofs", projector.get_n_dofs());
					table_H.add_value("ncells", projector.get_n_cells());
					table_H.add_value("L2", projector.get_L2_norm());
					table_H.add_value("H1", 0.0);
				}
				{ // Calculate B'.
					fname = dir + "solution_B_p" + std::to_string(p) +
						"_r" + std::to_string(r);

					table_B.add_value("r", r);
					table_B.add_value("p", p);

					if (SettingsSSOLIIIAXI::print_time_tables)
						std::cout << "Time table D \n";

					ExactSolutionSSOLIIIAXI_B exact_solution;

					ProjectAphiToBrz projector(
						p-1,
						problem.get_mapping_degree(),
						problem.get_tria(),
						problem.get_dof_handler(),
						problem.get_solution(),
						fname,
						& exact_solution,
						Settings::print_time_tables,
						Settings::project_exact_solution,
						Settings::log_cg_convergence);

					table_B.add_value("ndofs", projector.get_n_dofs());
					table_B.add_value("ncells", projector.get_n_cells());
					table_B.add_value("L2", projector.get_L2_norm()/mu_0);
					table_B.add_value("H1", 0.0);
				}
			}
			// Save convergence tables
			std::cout << "Table H\n";
			table_H.save(dir + "table_H_p" + std::to_string(p));

			std::cout << "Table B\n";
			table_B.save(dir + "table_B_p" + std::to_string(p));

		}
	}
};

int main()
{
	try
	{
		BatchSSOLIIIAXI batch;
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

