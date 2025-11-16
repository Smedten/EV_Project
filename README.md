# Building the project
### Disclaimer: This project has only been tried and tested on Ubuntu 18.04<= and Windows WSL.
A requirement for building the project is having the CPLEX solver installed on your PC. Instructions for this can be found here: https://www.leandro-coelho.com/install-cplex-on-linux-without-administration-privileges.

1. Git clone the repository.
    - Command: ```git clone https://github.com/Ravnholt7507/p9-projekt.git```.
2. Create a CPLEX solver environment variable.
    - Command: ```export CPLEX_STUDIO_DIR2211=/path/to/CPLEX_Studio_Community2211``` (For persistence, this line can also be added to your .bashrc).
3. Step into the build folder of the repository.
    - Command: ```cd build``` (assuming you are already in the project folder).
4. Create your own version of the makefile (as the one on GitHub might not be up to date).
    - Make sure CMake is installed on your computer, if not it can be installed with ```snap install cmake``` on Ubuntu systems.
    - Command: ```cmake ../CMakeLists.txt```.
5. Compile the project. 
    - Command: ```cmake --build .```.

# Running the project
1. To run the project, make sure the output program is executable.
    - Command: ```chmod +x ./output```.
2. Execute the program
    - Command: ```./output```

# Running unit tests
```./output --run_tests```

# Running the interactive demo
1. Ensure the Python dependencies are available
   - Command: ```pip install pandas matplotlib```
2. Launch the interactive dashboard from the `visuals` folder
   - Command: ```python3 visuals/demo_dashboard.py```
3. Use the radio buttons and slider to set aggregator and alignment modes, then press **Show demo** to update the plots.

# Generating graphs
To generate the graphs used in the report:
1. Step into the visuals folder (assuming you are in the projects root folder)
    - Command: ```cd visuals```
2. Make sure the pip libraries pandas and matplotlib are installed
    - Command: ```pip install pandas && pip install matplotlib```
3. Run the Python file
    - Command: ```python3 plots.py```
