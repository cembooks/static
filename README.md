<h1>Supplemented computer code</h1>

<h2> Introduction </h2>

This is the supplemental computer code for the book "Electromagnetics with
finite elements: static problems", ISBN 978-9090380933. The code contains
numerical experiments that illustrate application of the finite element
method to problems in electro- and magnetostatics. A number of boundary
value problems formulated in terms of electric scalar potential, total
magnetic scalar potential, reduced magnetic scalar potential, magnetic vector
potential, and current vector potential are solved numerically on two- and
three-dimensional problem domains.

The documentation of the computer code and the numerical experiments
exists in a form of online Logbook. The Logbook can be found at 
https://cembooks.nl.

<h2> Installation </h2>

To run the code do the following.

1) Get [deal.II](https://dealii.org) and install it.

2) Install texlive and GNU parallel.

3) Install a program for viewing vtk files. [Visit](https://visit.llnl.gov) 
software package of the Lawrence Livermore National Laboratory is a good option.

4) Change into the static/ and run

./clean  
./setup  
./build  
./run-all  

This will build and execute all numerical experiments. The last two steps will
take some time... Alternatively, you can build and run only the experiments you
are interested in. To do so, change into static/ and run

./clean  
./setup  

Then change into the directory of the experiment you are interested in, say
mms/. In this directory change into build/Release/ and run

./clean  
./build  

After that change into mms/bin/Release/ and run

./run-all  

This will execute the mms experiment. The results of the simulations are in
mms/bin/Release/Data/ directory.

