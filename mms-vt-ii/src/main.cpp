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
#include "project_Tz_to_Jxy.hpp"
#include "project_Axy_to_Bz.hpp"

#include "settings.hpp"

using namespace Misc;
using namespace StaticScalarSolver;
using namespace StaticVectorSolver;

/**
 * \brief This is a wrap-around class. It contains the main loop of the program
 * that implements the
 * [Method of manufactured solutions, vector potential (mms-vt-ii/)](@ref page_mms_vt_ii)
 * numerical experiment.
 *
 * The purpose of this class is to make the main function to be similar to the
 * main function of the deal.II
 * [Step-6](https://dealii.org/developer/doxygen/deal.II/step_6.html)
 * tutorial.
 *****************************************************************************/
class BatchMMSVTII : public SettingsMMSVTII
{
public:
	void run()
	{
		if (nr_threads_max > 0)
			MultithreadInfo::set_thread_limit(nr_threads_max);

#if DOMAIN__== 0
		std::string dir = "Data/circle/";
#endif

#if DOMAIN__ == 1
		std::string dir = "Data/square/";
#endif

		std::cout
			<< "Program: mss-vt-ii\n"
			<< "Dimensions: " << "2\n"
			<< "Domain type: " << DOMAIN__ << "\n"
			<< "Writing to: " << dir << "\n";

		MainOutputTable table_T(2);
		MainOutputTable table_J(2);
		MainOutputTable table_B(2);

		for (unsigned int p = 0; p < 3; p++)
		{
			table_T.clear();
			table_J.clear();
			table_B.clear();

			for (unsigned int r = 9; r < 13; r++)
			{
				table_T.add_value("r", r);
				table_T.add_value("p", p);

				table_J.add_value("r", r);
				table_J.add_value("p", p);

				table_B.add_value("r", r);
				table_B.add_value("p", p);

// Stage 0 -------------------------------------------------------------
				std::cout << "Stage 0: solving for T ...\n";

				SolverMMSVTII_T stage0(p+1, r,
					dir + "solution_T_p" + std::to_string(p) +
					"_r" + std::to_string(r));

				stage0.clear();
				table_T.add_value("ndofs", stage0.get_n_dofs());
				table_T.add_value("ncells", stage0.get_n_cells());
				table_T.add_value("L2", stage0.get_L2_norm());
				table_T.add_value("H1", stage0.get_H1_norm());

// Stage 1 -------------------------------------------------------------
				std::cout << "Stage 1: projecting T in H(curl) to Jf in H(div) ...\n";

				ExactSolutionMMSVTII_Jf exact_solution_Jf;

				ProjectTzToJxy<1> stage1(
					p,
					stage0.get_mapping_degree(),
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

				SolverMMSVTII_A stage2(
					p,
					stage0.get_tria(),
					stage0.get_dof_handler(),
					stage0.get_solution(),
				  r,
					dir + "solution_A_p" + std::to_string(p) +
					"_r" + std::to_string(r));

				stage2.clear();

// Stage 3 -------------------------------------------------------------
				std::cout << "Stage 3: projecting A in H(curl) to B in H(div) ...\n";

				ExactSolutionMMSVTII_B exact_solution_B;

				ProjectAxyToBz<3> stage3(
					p,
					1,
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
			std::cout << "Table T \n";
			table_T.save(dir + "main_table_T_p" + std::to_string(p));

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
		BatchMMSVTII batch;
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

