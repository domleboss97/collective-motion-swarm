#include <iostream>
#include <iomanip>
#include <cmath>    //cos and sin
#include <vector>   //used to track particles within a cell
#include <algorithm>//min and max
#include <iterator> //vector pointer operations
#include <ctime>    //time-based seed
#include <random>   //uniform random distribution
#include <fstream>

#define PI 3.14159265
#define MASS 1
#define part_diam 1.122
#define delta .0005
#define cutoff 2.5
#define cellsize 1.25
#define noise_strength 10
#define temp 10
#define Epsilon 1
#define strength 1

using namespace std;

//Contains data about each particle in simulation
class Swimmer{
    public:

		//Position
        double x_pos;
		double y_pos;

		//Sum of the forces
		double x_sum;
		double y_sum;

		Swimmer()
        {
            x_pos = 0;
            y_pos = 0;
            x_sum = 0;
            y_sum = 0;
        }

        void set_x_position(double x)
        {
            x_pos = x;
        }


        void set_y_position(double y)
        {
            y_pos = y;
        }

};


//Function which determines whether or not a particle is close enough to affect distance in vicsek model
double coord_distance(double x_1, double y_1, double x_2, double y_2)
{
    return sqrt(pow(x_1 - x_2, 2) + pow(y_1 - y_2, 2));
}

//unit vector calculations
double unit_x(double x_1, double y_1, double x_2, double y_2)
{
    return (x_1 - x_2)/ coord_distance(x_1, y_1, x_2, y_2);
}

double unit_y(double x_1, double y_1, double x_2, double y_2)
{
    return (y_1 - y_2)/ coord_distance(x_1, y_1, x_2, y_2);
}

//default_random_engine generator1(time(0)%91);
//default_random_engine generator2(time(0)%83);

int b, j, k, l, m, n, p; //For-loop counters
float t;
double x, y, dist, stage_x, stage_y;
double stage_pull_x, stage_pull_y, target_pull_x, target_pull_y;


int main()
{
    double time_track = 0;

    ofstream fout("positions.csv", ios::out);//| ios::app);
    ofstream fout3("Force_on_smart.csv", ios::out);

    int popsize, timesteps, fieldsize;
    double targetx_1, targety_1, targetx_2, targety_2, winteraction, interaction, angle;

    cout << "Population size?" << endl;
    cin >> popsize;

    cout << "Target 1 x?" << endl;
    cin >> targetx_1;

    cout << "Target 1 y?" << endl;
    cin >> targety_1;

    cout << "Target 2 x?" << endl;
    cin >> targetx_2;

    cout << "Target 2 y?" << endl;
    cin >> targety_2;

    cout << "Timesteps?" << endl;
    cin >> timesteps;

    fieldsize = 35;

    fout << popsize << endl << fieldsize << endl << temp << endl << targetx_1 << "," << targety_1 << endl << targetx_2 << "," << targety_2 << endl;
    fout << endl;

    fout << fixed;
    fout << setprecision(5);

    vector<Swimmer> parts(popsize);
    normal_distribution <double> pos_dance(0,1);

    //Initialization of particles
    for (b = 0; b < popsize - 1; b++)
    {
        double x = 1 + (floor(b/((fieldsize-2)/(cellsize)))) * cellsize;
        double y = 1 + fmod((b * cellsize), (fieldsize-2));

        parts[b].set_x_position(x);
        parts[b].set_y_position(y);

        //Prints the initial positions of the Particles
		fout << parts[b].x_pos << "," << parts[b].y_pos << endl;

    }

    b = popsize - 1;
    x = fieldsize - 4;
    y = .5 * fieldsize;

    parts[b].set_x_position(x);
    parts[b].set_y_position(y);

    fout << parts[b].x_pos << "," << parts[b].y_pos << endl;

	fout << endl;

	 for (t = 1; t <= timesteps; t += 1)
    {
        for (b = 0; b < popsize; b++)
        {
            //Loop though all the possible interactions
            for (l = 0; l < popsize; l++)
            {
                //If statement checks if particle is within 1 unit of another particle, and accounts for wraparound effect of board
                if (l != b  && (dist = coord_distance(parts[b].x_pos, parts[b].y_pos, parts[l].x_pos, parts[l].y_pos) < cutoff))
                {
                    //Interaction calculations
                    interaction = (12*Epsilon / part_diam) * ((pow(part_diam / dist, 13) - (pow(part_diam / dist, 7))));

                    // calculates the angle between the interacting particles
					angle = atan2(parts[l].y_pos - parts[b].y_pos, parts[l].x_pos - parts[b].x_pos);

                    //Sum of Interaction forces - separated into x and y components
					parts[b].x_sum += -1 * cos(angle) * (interaction);
					parts[j].x_sum += 1 * cos(angle) * (interaction);
                    parts[b].y_sum += -1 * sin(angle) * (interaction);
                    parts[j].y_sum += 1 * sin(angle) * (interaction);
                }
            }
            //Boundary Conditions
            //top
            if ((dist = coord_distance(parts[b].x_pos, parts[b].y_pos, parts[b].x_pos, fieldsize)) < cutoff)
            {
               winteraction = (12*Epsilon / part_diam) * (pow(part_diam / dist, 13));

               parts[b].y_sum += -1 * winteraction;
            }

            //bottom
            if ((dist = coord_distance(parts[b].x_pos, parts[b].y_pos, parts[b].x_pos, 0)) < cutoff)
            {
               winteraction = (12*Epsilon / part_diam) * (pow(part_diam / dist, 13));

               parts[b].y_sum += 1 * winteraction;
            }

            //right
            if ((dist = coord_distance(parts[b].x_pos, parts[b].y_pos, fieldsize, parts[b].y_pos)) < cutoff)
            {
               winteraction = (12*Epsilon / part_diam) * (pow(part_diam / dist, 13));

               parts[b].x_sum += -1 * winteraction;
            }

            //left
            if ((dist = coord_distance(parts[b].x_pos, parts[b].y_pos, 0, parts[b].y_pos)) < cutoff)
            {
               winteraction = (12*Epsilon / part_diam) * (pow(part_diam / dist, 13));

               parts[b].x_sum += 1 * winteraction;
            }
        }

        for (b = 0; b < popsize-1; b++)
        {
            stage_pull_x = strength * unit_x(parts[b].x_pos, parts[b].y_pos, stage_x, stage_y);
            stage_pull_y = strength * unit_y(parts[b].x_pos, parts[b].y_pos, stage_x, stage_y);

            parts[b].x_pos += (parts[b].x_sum + stage_pull_x) * delta;
            parts[b].y_pos += (parts[b].y_sum + stage_pull_y) * delta;

            // print statement
            fout << parts[b].x_pos << "," << parts[b].y_pos << endl;
        }

        b = popsize - 1;

        target_pull_x = strength * unit_x(parts[b].x_pos, parts[b].y_pos, stage_x, stage_y);
        target_pull_y = strength * unit_y(parts[b].x_pos, parts[b].y_pos, stage_x, stage_y);

        parts[b].x_pos += (parts[b].x_sum + target_pull_x) * delta;
        parts[b].y_pos += (parts[b].y_sum + target_pull_y) * delta;

        fout << parts[b].x_pos << "," << parts[b].y_pos << endl;

        fout << endl;

    }

    time_track = clock() - time_track;
    cout << "time = " << time_track/CLOCKS_PER_SEC << endl;

    return 0;
}
