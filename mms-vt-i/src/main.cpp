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

#include <deal.II/base/timer.h>
#include <deal.II/base/multithread_info.h>

#include <iostream>
#include <string>

#include "misc.hpp"
#include "solver.hpp"
#include "project_T_to_J.hpp"
#include "project_A_to_B.hpp"

#include "settings.hpp"

using namespace Misc;
using namespace StaticVectorSolver;

/**
 * \brief This is a wrap-around class. It contains the main loop of the program
 * that implements the
 * [Method of manufactured solutions, vector potential (mms-vt-i/)](@ref page_mms_vt_i)
 * numerical experiment.
 *
 * The purpose of this class is to make the main function to be similar to the
 * main function of the deal.II
 * [Step-6](https://dealii.org/developer/doxygen/deal.II/step_6.html)
 * tutorial.
 *****************************************************************************/
class BatchMMSVTI : public SettingsMMSVTI
{
public:
	void run()
	{
		if (nr_threads_max > 0)
			MultithreadInfo::set_thread_limit(nr_threads_max);


#if DOMAIN__ == 0
		const unsigned int mapping_degree = 1;
		const std::string dir = "Data/cube/";
#endif

#if DOMAIN__ == 1
		const unsigned int mapping_degree = 1;
		const std::string dir = "Data/sphere/";
#endif

		std::cout
			<< "Program: mss-vt-i\n"
			<< "Dimensions: " << "3\n"
			<< "Domain type: " << DOMAIN__ << "\n"
			<< "Writing to: " << dir << "\n";

		MainOutputTable table_J(3);
		MainOutputTable table_B(3);

		for (unsigned int p = 0; p < 3; p++)
		{
			table_J.clear();
			table_B.clear();

			for (unsigned int r = 5; r < 9; r++)
			{
				table_J.add_value("r", r);
				table_J.add_value("p", p);

				table_B.add_value("r", r);
				table_B.add_value("p", p);

// Stage 0 -------------------------------------------------------------
				std::cout << "Stage 0: solving for T ...\n";

				SolverMMSVTI_T stage0(p, mapping_degree, r,
					dir + "solution_T_p" + std::to_string(p) +
					"_r" + std::to_string(r));

				stage0.clear();

// Stage 1 -------------------------------------------------------------
				std::cout << "Stage 1: projecting T in H(curl) to Jf in H(div) ...\n";

				ExactSolutionMMSVTI_Jf exact_solution_Jf;

				ProjectTtoJ stage1(
					p,
					mapping_degree,
					stage0.get_tria(),
					stage0.get_dof_handler(),
					stage0.get_solution(),
					dir + "solution_J_p" + std::to_string(p) +
					"_r" + std::to_string(r),
					& exact_solution_Jf,
					Settings::print_time_tables,
					Settings::project_exact_solution,
					Settings::log_cg_convergence
					);

				table_J.add_value("ndofs", stage0.get_n_dofs());
				table_J.add_value("ncells", stage0.get_n_cells());
				table_J.add_value("L2", stage1.get_L2_norm());
				table_J.add_value("H1", 0.0);

				stage1.clear();

// Stage 2 --------------------------------------------------------------
				std::cout << "Stage 2: solving for A ...\n";

				SolverMMSVTI_A stage2(
					p,
					mapping_degree,
					stage0.get_tria(),
					stage0.get_dof_handler(),
					stage0.get_solution(),
				  r,
					dir + "solution_A_p" + std::to_string(p) +
					"_r" + std::to_string(r));

				stage2.clear();

// Stage 3 -------------------------------------------------------------
				std::cout << "Stage 3: projecting A in H(curl) to B in H(div) ...\n";

				ExactSolutionMMSVTI_B exact_solution_B;

				ProjectAtoB stage3(
					p,
					mapping_degree,
					stage0.get_tria(),
					stage2.get_dof_handler(),
					stage2.get_solution(),
					dir + "solution_B_p" + std::to_string(p) +
					"_r" + std::to_string(r),
					& exact_solution_B,
					Settings::print_time_tables,
					Settings::project_exact_solution,
					Settings::log_cg_convergence
					);

				table_B.add_value("ndofs", stage0.get_n_dofs());
				table_B.add_value("ncells", stage0.get_n_cells());
				table_B.add_value("L2", stage3.get_L2_norm());
				table_B.add_value("H1", 0.0);
			}

			std::cout << "Table J \n";
			table_J.save(dir + "main_table_J_p" + std::to_string(p));
			std::cout << "Table B \n";
			table_B.save(dir + "main_table_B_p" + std::to_string(p));
		}
	}
};

int main()
{
	try
	{
		BatchMMSVTI batch;
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

