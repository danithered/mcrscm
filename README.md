# Transition from mineral surfaces to vesicles at the origin of life

Origin of life spatio-temporal simulation of a 2D-3D transition from MCRS (surface metabolism) to SCM (cellularized metabolism). 
Simulation made in C++, plotting in R.
We are using ViennaRNA to map RNA sequences to secondary structures, to calculate dynamic properties. Replicators are tested in different ecological settings like associated to surfaces or in vesicules.

## Installation

### Prequisites

The simulations were compiled and run in Debian / Ubuntu environments.
In order to be able to compile the codes, the folowing softwares are needed:

- GNU Scientific Library
- Boost
- Vienna RNA package
- make, pkg-config, gcc

### Installation

- Download the code. It is suggested to download it with git:

        ```
        cd /your/path/
        git clone ...
        ```
- Compile the program: 
        
        ```
        cd /your/path/mcrscm
        make
        ```
        
        This will compile three programs: `mcrs` and `scm` for the corresponding spatial simulations, and `randseq` for the generation of random sequence pool

## Usage

- Choose a ruleset for the enzymatic activities. 
    - The default ones, that have been used for the article, is set for deafult. If it is your choice, you do not have to do anything, jump to the next point.
    - If you want to use an other set, drawn from our way of generating method, do the following:   
        1) `make gen`
        2) run `./gen` to overwrite the original with an alternative mapping in "IN/str" or `./gen custom/path/` to put into "custom/path/"
        3) run `src/shuffle.sh IN/str/` to shuffle the order of rules
        4) run `src/create_mapping.sh N > IN/str/mapping_AN.txt`, where `N` stands for the number of enzymatic activities in a simulation.
    - If you want to define rules for yourself, edit a *IN/str/mappingAN.txt* file, according to the followings: 
        - first  line sets the number of enzymatic activities or system size $A$
        - after the first line, there are blocks of motif definitions
        - in each bock the first line describes a spatial pattern to look for. The spatial pattern is in a dot-bracket format. Before the spetial pattern definition there is a number defining how many alternate sequenc patterns to look for.
        - after the spatial pattern the sequence patterns are described intwo lines definitions. The first line sets the sequences in  a format of sequences of position - base pair identity pairs. The second line renders activities to the enzymatic activties. E.g. in case of $A=3$, this line looks like "$\alpha_1$ $\alpha_2$ $\alpha_3$", where each $\alpha$ corresponds to the first, second an third enzymatic activities. In the original article in the motif catalysed enzymatic activity 3, we used "0 0 1".
        - Please use maximum 3 long sequence definitions. In this case the following transformation will be applied to define activity $a(n)$ for motifs, where $n$ is the number of sequnce rules fulfilled: $a(0)=0.0$, $a(1)=0.1 \alpha$, $a(2)=0.8 \alpha$ and $a(3) = \alpha$.
- Generate a pool of random replicators. By default the program generates 10 000 000 sequences in a length of $\lambda=45$ in a Poisson distribution. If you wish to change that please refer to *src/randseq.cpp*! To start generation type: `./randseq --par_str_pool IN/str/mappingA10.txt`. Please note, that this step will take a while!
- Run the simulations. `./mcrs` and `./scm` starts the corresponding simulations. To set the parameters, please refer to **Model parameters** section!

## Model parameters

Parameters can be included in any order after the program name (e.g. `./mcrs --ncol 80 --nrow 90` and `./mcrs --nrow 90 --ncol 80`  will refer to the same parameter settings.) 

### MCRS parameters

|  Parameter              | Description                                                                                                       | Default value          |
|-------------------------|-------------------------------------------------------------------------------------------------------------------|------------------------|
| `--par_maxtime`         | The maximum duration of a simulation                                                                              | 2000000                |
| `--par_ncol`            | Number of columns of the grid                                                                                     |     300                |
| `--par_nrow`            | Number of rows of the grid                                                                                        |     300                | 
| `--par_seed`            | Random seed. If -1 than the system will be the seed.                                                              |    -1                  |
| `--par_output_interval` | Output interval.                                                                                                  | 1000                   |
| `--par_save_interval`   | The interval when the full grid is saved.                                                                         | 1000000                |
| `--par_ID`              | The ID of the simualtion. The forlder name where the results will be saved.                                       | "test\0"               |
| `--par_str_pool`        | The mapping of enzymatic activities. It willset the system size too.                                              | "IN/str/mappingA3.txt" |
| `--par_load`            | Loading simulation from a grid.                                                                                   |                        | 
| `--par_diffusion_rate`  | Diffusion rate                                                                                                    | 4                      |
| `--par_claimEmpty`      | Claim of empty place to remain empty                                                                              | 0.1                    |
| `--par_Nmet`            | Size of metabolic neighbourhood. Determined in square of euclidean distance. For vonNeumann use 1, for Moore 2.   | 4                      | 
| `--par_Nrep`            | Size of replication neighbourhood. Determined in square of euclidean distance. For vonNeumann use 1, for Moore 2. | 4                      |
| `--par_substitution`    | Per base probability of substituion mutations.                                                                    | 0.005                  |
| `--par_insertion`       | Per base probability of insertion mutations.                                                                      | 0.0005                 |
| `--par_deletion`        | Per base probability of deletion mutations.                                                                       | 0.0005                 |
| `--par_c`               | Minus of $c$ parameter. See article Eq. 2.                                                                        | -0.3                   |
| `--par_sigma`           | $\sigma$ parameter. See article Eq. 5.                                                                            | 1.1                    |
| `--par_gc_bonus`        | Parameter setting additional GC bases’ contribution to enzymatic activities.                                      | -0.3                   |
| `--par_g`               | $g$ parameter. See article Eq. 3.                                                                                 | 10                     |
| `--par_b1`              | $b_1$ parameter. See article Eq. 3.                                                                               | 0.75                   |
| `--par_b2`              | $b_2$ parameter. See article Eq. 3.                                                                               | 0.005                  |
| `--par_ll`              | Parameter $l$ +1. See article Eq. 3.                                                                              | 2                      | 
| `--par_Emin`            | $E_{min}$ - Minimal folding energy                                                                                | -25.0                  |
| `--par_flexPdeg`        | $d$ parameter. See article Eq. 1.                                                                                 | 0.2                    |
| `--par_bubble_interval` | Frequency of sampling vesicles from the grid.                                                                     | 0                      |
| `--par_no_bubi`         | Number of vesicles sampled from the grid per generation.                                                          | 0                      |
| `--par_mean_bubblesize` | Average vesicle size.                                                                                             | 6.0                    |
| `--par_sd_bubblesize`   | SD of vesicle size.                                                                                               | 1                      |

### SCM parameters

|  Parameter                | Description                                                                    | Default value          |
|---------------------------|--------------------------------------------------------------------------------|------------------------|
| `--par_maxtime`           | The maximum duration of a simulation                                           | 2000000                |
| `--par_num_input_content` | Number of replicators in initial vesicles.                                     |   25                   |
| `--par_splitfrom`         | Splitsize ($S$)                                                                |   50                   |
| `--par_poolsize`          | Number of vesicles ($N$)                                                       |   300                  |
| `--par_seed`              | Random seed. If -1 than the system will be the seed.                           |    -1                  |
| `--par_output_interval`   | Output interval.                                                               | 1000                   |
| `--par_save_interval`     | The interval when the full grid is saved.                                      | 1000000                |
| `--par_ID`                | The ID of the simualtion. The forlder name where the results will be saved.    | "test\0"               |
| `--par_str_pool`          | The mapping of enzymatic activities. It willset the system size too.           | "IN/str/mappingA3.txt" |
| `--par_load`              | Loading simulation from a grid.                                                |                        | 
| `--par_bubbles`           | Location (folder) of vesicles to initialise simulation from.                   | "\0"                   | 
| `--par_substitution`      | Per base probability of substituion mutations.                                 | 0.005                  |
| `--par_insertion`         | Per base probability of insertion mutations.                                   | 0.0005                 |
| `--par_deletion`          | Per base probability of deletion mutations.                                    | 0.0005                 |
| `--par_c`                 | Minus of $c$ parameter. See article Eq. 2.                                     | -0.3                   |
| `--par_sigma`             | $\sigma$ parameter. See article Eq. 5.                                         | 1.1                    |
| `--par_gc_bonus`          | Parameter setting additional GC bases’ contribution to enzymatic activities.   | -0.3                   |
| `--par_g`                 | $g$ parameter. See article Eq. 3.                                              | 10                     |
| `--par_b1`                | $b_1$ parameter. See article Eq. 3.                                            | 0.75                   |
| `--par_b2`                | $b_2$ parameter. See article Eq. 3.                                            | 0.005                  |
| `--par_ll`                | Parameter $l$ +1. See article Eq. 3.                                           | 2                      | 
| `--par_Emin`              | $E_{min}$ - Minimal folding energy                                             | -25.0                  |
| `--par_flexPdeg`          | $d$ parameter. See article Eq. 1.                                              | 0.2                    |

## Results

We have shared the results of Fig. 1. in the *results* folder. As the all our results spans to terrabites of files, we choose to share only the most substanitial results. Other results can be regenerated by our codes (please note that these are stochastic simulations!)


