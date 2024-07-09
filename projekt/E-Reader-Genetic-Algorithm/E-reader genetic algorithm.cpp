# define _CRT_SECURE_NO_WARNINGS

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>
# include <cstring>
# include <vector>
# include "CImg.h"

using namespace cimg_library;
using namespace std;
// 
//  Change any of these parameters to match your needs 
//
# define MAXGENS 70
# define NVARS 8
# define PMUTATION 0.05
# define POPSIZE 200
# define PXOVER 0.8
//
//  Each GENOTYPE is a member of the population, with
//  gene: a string of variables,
//  fitness: the fitness value
//  upper: the variable upper bounds,
//  lower: the variable lower bounds,
//  rfitness: the relative fitness,
//  cfitness: the cumulative fitness
//
struct genotype {
    double gene[NVARS];
    double fitness;
    double upper[NVARS];
    double lower[NVARS];
    double rfitness;
    double cfitness;
};

struct genotype population[POPSIZE + 1];
struct genotype newpopulation[POPSIZE + 1];
//
//  All program functions listed here
//
int main();
void crossover(int& seed);
void elitist();
void evaluate(vector<int> data);
int i4_uniform_ab(int a, int b, int& seed);
void initialize(int& seed);
void keep_the_best();
void mutate(int& seed);
double r8_uniform_ab(double a, double b, int& seed);
void report(int generation);
void selector(int& seed);
void timestamp();
void Xover(int one, int two, int& seed);
bool correct_gene(double* gene);

//****************************************************************************80

int main() {

//****************************************************************************80
//
//  Purpose:
//
//    MAIN supervises the genetic algorithm.
//
//  Discussion:
//
//    Each generation involves selecting the best 
//    members, performing crossover & mutation and then 
//    evaluating the resulting population, until the terminating 
//    condition is satisfied   
//
//    This is a more complex genetic algorithm implementation, 
//    where the evaluation function is minimized and  
//    fitness of an individual is the same as the value
//    of the mentioned evaluation function.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//    17 May 2024
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//    This C++ updated version by Wojciech Michaluk.
//
//  Reference:
//
//    Zbigniew Michalewicz,
//    Genetic Algorithms + Data Structures = Evolution Programs,
//    Third Edition,
//    Springer, 1996,
//    ISBN: 3-540-60676-9,
//    LC: QA76.618.M53.
//
//  Parameters:
//
//    MAXGENS is the maximum number of generations.
//    NVARS is the number of problem variables.
//    PMUTATION is the probability of mutation.
//    POPSIZE is the population size. 
//    PXOVER is the probability of crossover.                          
//
    CImg<int> image;
    string filepath = "images/img1.bmp";
    //
    //  For other example files:
    //  string filepath = "images/img2.bmp";
    //  string filepath = "images/img3.bmp";
    //  string filepath = "images/img4.bmp";
    //
    int generation;
    int i;
    int seed;

    image.load(filepath.c_str());
    int *data = image.data();
    vector<int> pixels_data;
    vector<int> original_pixels_data;

    for (int i = 0; i < image.width() * image.height(); i++) {
        pixels_data.push_back(data[i]);
        original_pixels_data.push_back(data[i]);
    }

    sort(pixels_data.begin(), pixels_data.end());

    timestamp();
    cout << "\n";
    cout << "E-reader genetic algorithm:\n";
    cout << "  C++ version\n";
    cout << "  An example of a genetic algorithm in action.\n";

    if (NVARS < 2) {
        cout << "\n";
        cout << "  The crossover modification will not be available,\n";
        cout << "  since it requires 2 <= NVARS.\n";
    }

    seed = 123456789;

    initialize(seed);
    evaluate(pixels_data);
    keep_the_best();

    for (generation = 0; generation < MAXGENS; generation++) {
        selector(seed);
        crossover(seed);
        mutate(seed);
        report(generation);
        evaluate(pixels_data);
        elitist();
    }

    cout << "\n";
    cout << "  Best member after " << MAXGENS << " generations:\n";
    cout << "\n";

    for (i = 0; i < NVARS; i++) {
        cout << "  var(" << i << ") = " << population[POPSIZE].gene[i] << "\n";
    }

    cout << "\n";
    cout << "  Best fitness = " << population[POPSIZE].fitness << "\n";
    //
    //  Display generated image and then terminate
    //
    CImg<int> new_image(image.width(), image.height());
    vector<int> brightness;
    double best_gene[NVARS];
    double assigned_new_values[256];

    int index = 0;
    int pixel_value = pixels_data.at(0);
    int color;

    for (int i = 0; i < NVARS; i++) {
        best_gene[i] = population[POPSIZE].gene[i];
    }

    for (int i = 0; i < NVARS - 1; i++) {
        while (index < pixels_data.size() && abs(pixel_value - best_gene[i]) < abs(pixel_value - best_gene[i + 1])) {
            assigned_new_values[pixel_value] = best_gene[i];
            index++;

            if (index != pixels_data.size()) {
                pixel_value = pixels_data.at(index);
            }
        }
    }

    while (index < pixels_data.size()) {
        assigned_new_values[pixel_value] = best_gene[NVARS - 1];
        index++;

        if (index != pixels_data.size()) {
            pixel_value = pixels_data.at(index);
        }
    }

    for (int i = 0; i < image.width(); i++) {
        for (int j = 0; j < image.height(); j++) {
            color = (int)assigned_new_values[original_pixels_data[j * image.width() + i]];
            new_image.fillC(i, j, 0, color, color, color);
        }
    }

    new_image.display();
    //
    //  In order to save an image:
    //  new_image.save_bmp("new_img<image_number>.bmp");
    //
    cout << "\n";
    cout << "E-reader genetic algorithm:\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp();

    return 0;
}
//****************************************************************************80

void crossover(int& seed) {

//****************************************************************************80
// 
//  Purpose:
//
//    CROSSOVER selects two parents for the single point crossover.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, int FIRST, is a count of the number of members chosen.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
    const double a = 0.0;
    const double b = 1.0;
    int mem;
    int one;
    int first = 0;
    double x;

    for (mem = 0; mem < POPSIZE; mem++) {
        x = r8_uniform_ab(a, b, seed);

        if (x < PXOVER) {
            if (++first % 2 == 0) {
                Xover(one, mem, seed);
            } else {
                one = mem;
            }
        }
    }
}
//****************************************************************************80

void elitist() {

//****************************************************************************80
// 
//  Purpose:
//
//    ELITIST stores the best member of the previous generation.
//
//  Discussion:
//
//    The best member of the previous generation is stored as 
//    the last in the array. If the best member of the current 
//    generation is worse then the best member of the previous 
//    generation, the latter one would replace the worst member 
//    of the current population.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//    17 May 2024
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//    This C++ updated version by Wojciech Michaluk.
//
//  Local parameters:
//
//    Local, double BEST, the best fitness value.
//    Local, double WORST, the worst fitness value.
//
    int i;
    double best = population[0].fitness;
    int best_mem;
    double worst = population[0].fitness;
    int worst_mem;

    for (i = 0; i < POPSIZE - 1; i++) {
        if (population[i + 1].fitness > population[i].fitness) {
            if (best >= population[i].fitness) {
                best = population[i].fitness;
                best_mem = i;
            }

            if (population[i + 1].fitness >= worst) {
                worst = population[i + 1].fitness;
                worst_mem = i + 1;
            }
        } else {
            if (population[i].fitness >= worst) {
                worst = population[i].fitness;
                worst_mem = i;
            }

            if (best >= population[i + 1].fitness) {
                best = population[i + 1].fitness;
                best_mem = i + 1;
            }
        }
    }
    // 
    //  If the best individual from the new population is better than 
    //  the best individual from the previous population, then 
    //  copy the best from the new population; else replace the 
    //  worst individual from the current population with the 
    //  best one from the previous generation                     
    //
    if (population[POPSIZE].fitness >= best) {
        for (i = 0; i < NVARS; i++) {
            population[POPSIZE].gene[i] = population[best_mem].gene[i];
        }

        population[POPSIZE].fitness = population[best_mem].fitness;
    } else {
        for (i = 0; i < NVARS; i++) {
            population[worst_mem].gene[i] = population[POPSIZE].gene[i];
        }

        population[worst_mem].fitness = population[POPSIZE].fitness;
    }
}
//****************************************************************************80

void evaluate(vector<int> data) {

//****************************************************************************80
// 
//  Purpose:
//
//    EVALUATE implements the user-defined evaluation function.
//
//  Discussion:
//
//    Each time this is changed, the code has to be recompiled.
//    The current function is calculating the sum of squares of
//    distances (differences) between original image pixels 
//    brightness and those of current population members.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//    17 May 2024
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//    This C++ updated version by Wojciech Michaluk.
//
//  Parameters:
//    Input, vector<int> DATA, containing information about proposed pixels brightness.
//
    int member;
    int i;
    int index;
    int pixel_value;
    double current_gene;
    double next_gene;
    double squares_sum;

    for (member = 0; member < POPSIZE; member++) {
        squares_sum = 0;
        index = 0;
        pixel_value = data.at(0);

        for (i = 0; i < NVARS - 1; i++) {
            current_gene = population[member].gene[i];
            next_gene = population[member].gene[i + 1];

            while (index < data.size() && abs(pixel_value - current_gene) < abs(pixel_value - next_gene)) {
                squares_sum += (pixel_value - current_gene) * (pixel_value - current_gene);
                index++;

                if (index != data.size()) {
                    pixel_value = data.at(index);
                }
            }
        }

        while (index < data.size()) {
            squares_sum += (pixel_value - population[member].gene[NVARS - 1]) * (pixel_value - population[member].gene[NVARS - 1]);
            index++;

            if (index != data.size()) {
                pixel_value = data.at(index);
            }
        }

        population[member].fitness = squares_sum;
    }
}
//****************************************************************************80

int i4_uniform_ab(int a, int b, int& seed) {

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int A, B, the limits of the interval.
//    Input/output, int &SEED, the "seed" value, which should NOT be 0. On output, SEED has been updated.
//    Output, int I4_UNIFORM, a number between A and B.
//
    int c;
    const int i4_huge = 2147483647;
    int k;
    float r;
    int value;

    if (seed == 0) {
        cerr << "\n";
        cerr << "I4_UNIFORM_AB - Fatal error!\n";
        cerr << "  Input value of SEED = 0.\n";
        exit(1);
    }
    //
    //  Guarantee A <= B
    //
    if (b < a) {
        c = a;
        a = b;
        b = c;
    }

    k = seed / 127773;
    seed = 16807 * (seed - k * 127773) - k * 2836;

    if (seed < 0) {
        seed = seed + i4_huge;
    }

    r = (float)(seed) * 4.656612875E-10;
    //
    //  Scale R to lie between A - 0.5 and B + 0.5
    //
    r = (1.0 - r) * ((float)a - 0.5) + r * ((float)b + 0.5);
    //
    //  Use rounding to convert R to an integer between A and B
    //
    value = round(r);
    //
    //  Guarantee A <= VALUE <= B
    //
    if (value < a) {
        value = a;
    } else if (b < value) {
        value = b;
    }

    return value;
}
//****************************************************************************80

void initialize(int& seed) {

//****************************************************************************80
// 
//  Purpose:
//
//    INITIALIZE initializes the genes within the variables bounds. 
//
//  Discussion:
//
//    It also initializes (to DBL_MAX) all fitness values for each
//    member of the population. It randomly generates values between
//    bounds of 0 (inclusive) and 256 (exclusive) for each gene of 
//    each genotype in the population.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//    17 May 2024
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//    This C++ updated version by Wojciech Michaluk.
//
//  Parameters:
//
//    Input, string FILENAME, the name of the input file.
//    Input/output, int &SEED, a seed for the random number generator.
//
    int lbound = 0;
    int ubound = 256;
    vector<int> all_values;
    vector<int> vars;

    for (int i = lbound; i < ubound; i++) {
        all_values.push_back(i);
    }

    for (int i = 0; i < POPSIZE; i++) {
        population[i].fitness = 0;
        population[i].rfitness = 0;
        population[i].cfitness = 0;

        random_shuffle(all_values.begin(), all_values.end());
        vars.push_back(lbound);

        for (int j = 0; j < NVARS; j++) {
            vars.push_back(all_values.at(j));
        }

        sort(vars.begin(), vars.end());
        vars.push_back(ubound);

        population[i].lower[0] = lbound;
        population[i].upper[NVARS - 1] = ubound;

        for (int j = 0; j < NVARS; j++) {
            population[i].gene[j] = vars.at(j);
        }

        vars.clear();
    }

    population[POPSIZE].fitness = DBL_MAX;
}
//****************************************************************************80

void keep_the_best() {

//****************************************************************************80
// 
//  Purpose:
//
//    KEEP_THE_BEST keeps track of the best member of the population. 
//
//  Discussion:
//
//    Note that the last entry in the array population holds a 
//    copy of the best individual.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, int CUR_BEST, the index of the best individual.
//
    int cur_best = 0;
    int mem;
    int i;

    for (mem = 0; mem < POPSIZE; mem++) {
        if (population[POPSIZE].fitness > population[mem].fitness) {
            cur_best = mem;
            population[POPSIZE].fitness = population[mem].fitness;
        }
    }
    // 
    //  Once the best member in the population is found, copy the genes
    //
    for (i = 0; i < NVARS; i++) {
        population[POPSIZE].gene[i] = population[cur_best].gene[i];
    }
}
//****************************************************************************80

void mutate(int& seed) {

//****************************************************************************80
// 
//  Purpose:
//
//    MUTATE performs a random uniform mutation. 
//
//  Discussion:
//
//    A variable selected for mutation is replaced by a random value 
//    between the lower and upper bounds of this variable, which for
//    genes "in the middle" (all but far left and far right) is
//    equivalent to values between previous gene and next gene, as
//    after the mutation happens, the order should be kept.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//    17 May 2024
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//    This C++ updated version by Wojciech Michaluk.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
    const double a = 0.0;
    const double b = 1.0;
    int i;
    int j;
    double lbound;
    double ubound;
    double x;

    for (i = 0; i < POPSIZE; i++) {
        for (j = 0; j < NVARS; j++) {
            x = r8_uniform_ab(a, b, seed);

            if (x < PMUTATION) {
                //
                //  Assuming NVARS >= 2
                //
                if (j == 0) {
                    lbound = population[i].lower[j];
                    ubound = population[i].gene[j + 1];
                    population[i].gene[j] = i4_uniform_ab(lbound, ubound - 1, seed);
                } else if (j == NVARS - 1) {
                    lbound = population[i].gene[j - 1];
                    ubound = population[i].upper[j];
                    population[i].gene[j] = i4_uniform_ab(lbound + 1, ubound - 1, seed);
                } else {
                    lbound = population[i].gene[j - 1];
                    ubound = population[i].gene[j + 1];
                    population[i].gene[j] = i4_uniform_ab(lbound + 1, ubound - 1, seed);
                }
            }
        }
    }
}
//****************************************************************************80

double r8_uniform_ab(double a, double b, int& seed) {

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_AB returns a scaled pseudorandom R8.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the limits of the interval.
//    Input/output, int &SEED, the "seed" value, which should NOT be 0. On output, SEED has been updated.
//    Output, double R8_UNIFORM_AB, a number strictly between A and B.
//
    int i4_huge = 2147483647;
    int k;
    double value;

    if (seed == 0) {
        cerr << "\n";
        cerr << "R8_UNIFORM_AB - Fatal error!\n";
        cerr << "  Input value of SEED = 0.\n";
        exit(1);
    }

    k = seed / 127773;
    seed = 16807 * (seed - k * 127773) - k * 2836;

    if (seed < 0) {
        seed = seed + i4_huge;
    }

    value = (double)(seed) * 4.656612875E-10;
    value = a + (b - a) * value;

    return value;
}
//****************************************************************************80

void report(int generation) {

//****************************************************************************80
// 
//  Purpose:
//
//    REPORT reports progress of the simulation. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, double avg, the average population fitness.
//    Local, double best_val, the best population fitness.
//    Local, double square_sum, square of sum for std calc.
//    Local, double stddev, standard deviation of population fitness.
//    Local, double sum, the total population fitness.
//    Local, double sum_square, sum of squares for std calc.
//
    double avg;
    double best_val;
    int i;
    double square_sum;
    double stddev;
    double sum;
    double sum_square;

    if (generation == 0) {
        cout << "\n";
        cout << "  Generation       Best            Average       Standard \n";
        cout << "  number           value           fitness       deviation \n";
        cout << "\n";
    }

    sum = 0.0;
    sum_square = 0.0;

    for (i = 0; i < POPSIZE; i++) {
        sum += population[i].fitness;
        sum_square += population[i].fitness * population[i].fitness;
    }

    avg = sum / (double)POPSIZE;
    square_sum = avg * avg * POPSIZE;
    stddev = sqrt((sum_square - square_sum) / (POPSIZE - 1));
    best_val = population[POPSIZE].fitness;

    cout << "  " << setw(8) << generation
         << "  " << setw(14) << best_val
         << "  " << setw(14) << avg
         << "  " << setw(14) << stddev << "\n";
}
//****************************************************************************80

void selector(int& seed) {

//****************************************************************************80
// 
//  Purpose:
//
//    SELECTOR is the selection function.
//
//  Discussion:
//
//    Standard proportional selection for minimization problems incorporating 
//    the elitist model. This makes sure that the best member always survives.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//    17 May 2024
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//    This C++ updated version by Wojciech Michaluk.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
    const double a = 0.0;
    const double b = 1.0;
    int i;
    int j;
    int mem;
    double p;
    double sum = 0.0;
    //
    //  Find the total fitness of the population
    //
    for (mem = 0; mem < POPSIZE; mem++) {
        sum += population[mem].fitness;
    }
    //
    //  Calculate the relative fitness of each member
    //
    for (mem = 0; mem < POPSIZE; mem++) {
        population[mem].rfitness = (sum - population[mem].fitness) / (sum * (double) (POPSIZE - 1));
    }
    // 
    //  Calculate the cumulative fitness
    //
    population[0].cfitness = population[0].rfitness;
    
    for (mem = 1; mem < POPSIZE; mem++) {
        population[mem].cfitness = population[mem - 1].cfitness + population[mem].rfitness;
    }
    // 
    //  Select survivors using cumulative fitness
    //
    for (i = 0; i < POPSIZE; i++) {
        p = r8_uniform_ab(a, b, seed);

        if (p < population[0].cfitness) {
            newpopulation[i] = population[0];
        } else {
            for (j = 0; j < POPSIZE; j++) {
                if (population[j].cfitness <= p && p < population[j + 1].cfitness) {
                    newpopulation[i] = population[j + 1];
                }
            }
        }
    }
    // 
    //  Overwrite the old population with the new one
    //
    for (i = 0; i < POPSIZE; i++) {
        population[i] = newpopulation[i];
    }
}
//****************************************************************************80

void timestamp() {

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 October 2003
//    17 May 2024
//
//  Author:
//
//    John Burkardt
//    Updated C++ version by Wojciech Michaluk.
//
    const int TIME_SIZE = 40;
    static char time_buffer[TIME_SIZE];
    time_t now = time(NULL);
    const struct tm* tm = localtime(&now);

    strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);
    cout << time_buffer << "\n";
}
//****************************************************************************80

void Xover(int one, int two, int& seed) {

//****************************************************************************80
// 
//  Purpose:
//
//    XOVER performs crossover of the two selected parents. If genes are not
//    ordered after the crossover (which is required), then they are swapped.
//    If the order still is not maintained, the crossover is reversed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//    17 May 2024
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//    This C++ updated version by Wojciech Michaluk.
//
//  Local parameters:
//
//    Local, int point, the crossover point.
//    Local, double[] mem_one and mem_two storing parents' genes
//    so that they can be restored if crossover is reversed.
//
//  Parameters:
//
//    Input, int ONE, TWO, the indices of the two parents.
//    Input/output, int &SEED, a seed for the random number generator.
//
    int i;
    double t;
    double mem_one [NVARS];
    double mem_two [NVARS];
    // 
    //  Select the crossover point
    //
    int point = i4_uniform_ab(0, NVARS - 1, seed);
    //
    //  Save current genes if reversing is needed
    //
    for (i = 0; i < NVARS; i++) {
        mem_one[i] = population[one].gene[i];
        mem_two[i] = population[two].gene[i];
    }
    //
    //  Swap genes in positions 0 through POINT-1
    //
    for (i = 0; i < point; i++) {
        t = population[one].gene[i];
        population[one].gene[i] = population[two].gene[i];
        population[two].gene[i] = t;
    }
    //
    //  Reverse if operation is `illegal`, i. e. assigned values are not ordered
    //  even after swapping genes at the crossover point
    //
    if (point > 0) {
        if (population[one].gene[point - 1] > population[one].gene[point]) {
            t = population[one].gene[point - 1];
            population[one].gene[point - 1] = population[one].gene[point];
            population[one].gene[point] = t;

            if (!correct_gene(population[one].gene)) {
                for (int i = 0; i < NVARS; i++) {
                    population[one].gene[i] = mem_one[i];
                    population[two].gene[i] = mem_two[i];
                }

                return;
            }
        }

        if (population[two].gene[point - 1] > population[two].gene[point]) {
            t = population[two].gene[point - 1];
            population[two].gene[point - 1] = population[two].gene[point];
            population[two].gene[point] = t;

            if (!correct_gene(population[two].gene)) {
                for (int i = 0; i < NVARS; i++) {
                    population[one].gene[i] = mem_one[i];
                    population[two].gene[i] = mem_two[i];
                }

                return;
            }
        }
    } 
    
    if (point > 1) {
        if (population[one].gene[point - 1] < population[one].gene[point - 2]) {
            t = population[one].gene[point - 1];
            population[one].gene[point - 1] = population[one].gene[point - 2];
            population[one].gene[point - 2] = t;

            if (!correct_gene(population[one].gene)) {
                for (int i = 0; i < NVARS; i++) {
                    population[one].gene[i] = mem_one[i];
                    population[two].gene[i] = mem_two[i];
                }

                return;
            }
        }

        if (population[two].gene[point - 1] < population[two].gene[point - 2]) {
            t = population[two].gene[point - 1];
            population[two].gene[point - 1] = population[two].gene[point - 2];
            population[two].gene[point - 2] = t;

            if (!correct_gene(population[two].gene)) {
                for (int i = 0; i < NVARS; i++) {
                    population[one].gene[i] = mem_one[i];
                    population[two].gene[i] = mem_two[i];
                }

                return;
            }
        }
    }
}
//****************************************************************************80

bool correct_gene(double* gene) {

//****************************************************************************80
// 
//  Purpose:
//
//    CORRECT_GENE checks if the gene is correct, which means it
//    checks if all values are ordered ascending, without repetitions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2024
//
//  Author:
//
//    This C++ version by Wojciech Michaluk.
//
//  Parameters:
//
//    Input, double* GENE, the gene to be checked if correct.
//
    double prev = gene[0];

    for (int i = 1; i < NVARS; i++) {
        if (gene[i] <= prev) {
            return false;
        }

        prev = gene[i];
    }

    return true;
}
//****************************************************************************80
