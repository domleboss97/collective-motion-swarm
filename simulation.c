#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define TO_RADIANS(degrees)((M_PI / 180.0) * degrees);

typedef struct agent agent;
typedef struct square square;

struct agent {
  int agent_num;
  double x_coord;
  double y_coord;
  double x_dir;
  double y_dir;
};

struct square {
  int *active_agents;
  int index;
};

int in_radius(agent *first_agent, agent *second_agent);
double calc_dist(double x_1, double y_1, double x_2, double y_2);
 
//Initialize global variables
double NOISE = 0.10;
int BOARD_LENGTH = 100;
int AGENT_NUMBER = 10000;
double ORDER_PARAM = 0.0;
int COUNTER = 0;
square **BOARD;

//Calculate the order parameter based on the orientation of the agents
double order_param(agent *agent_list) {
  double avg_x = 0.0;
  double avg_y = 0.0;
  for (int i = 0; i < AGENT_NUMBER; i++) {
    avg_x = avg_x + agent_list[i].x_dir;
    avg_y = avg_y + agent_list[i].y_dir;
  }
  double total = sqrt(pow(avg_x/(double) AGENT_NUMBER, 2) + pow(avg_y / (double) AGENT_NUMBER, 2));
  //return (total / (double) AGENT_NUMBER);
  return total;
}

//Update the orientaton of the agents to align with its neighbors
void update_orientation(agent *next_agent, agent *agent_list) {
  int x_cell = (int) floor(next_agent->x_coord);
  int y_cell = (int) floor(next_agent->y_coord);
  int neighbor_count = 0;
  double new_x = 0.0;
  double new_y = 0.0;

  //Loop through only surrounding squares to optimize search time
  for (int x = -1; x <= 1; x++) {
    for (int y = -1; y <= 1; y++) {
      int x_tmp = x_cell + x;
      int y_tmp = y_cell + y;
      if (x_tmp < 0) {
	       x_tmp = x_tmp + BOARD_LENGTH;
      } 
      if (y_tmp < 0) {
	       y_tmp = y_tmp + BOARD_LENGTH;
      }
      if (x_tmp >= BOARD_LENGTH) {
	       x_tmp = x_tmp - BOARD_LENGTH;
      }
      if (y_tmp >= BOARD_LENGTH) {
	       y_tmp = y_tmp - BOARD_LENGTH;
      }
      //Loop through agents on given square
      for (int i = 0; i <= BOARD[x_tmp][y_tmp].index; i++) {
	       if (!(BOARD[x_tmp][y_tmp].active_agents[i] == next_agent->agent_num)) {
	         int neighbor_number = BOARD[x_tmp][y_tmp].active_agents[i];
	         if (in_radius(next_agent, &agent_list[neighbor_number])) {
	           neighbor_count++;
	           new_x = new_x + agent_list[neighbor_number].x_dir;
	           new_y = new_y + agent_list[neighbor_number].y_dir;
	         }
	       } else {
         }
      }

    }
  }
  
  //Incorporate angular noise into final direction
  double noise_angle = (((double) rand() / (double) RAND_MAX) * 360) - 180;
  noise_angle = TO_RADIANS(noise_angle);
  double current_angle;

  if (neighbor_count > 0) {
    current_angle = atan2(new_y, new_x);
  } else {
    current_angle = atan2(next_agent->y_dir, next_agent->x_dir);
  }
  current_angle = (NOISE * noise_angle) + current_angle;
  next_agent->x_dir = cos(current_angle);
  next_agent->y_dir = sin(current_angle);
  
}

//Check whether two particles are within each other's radii
int in_radius(agent *first_agent, agent *second_agent) {
  double x_1 = first_agent->x_coord;
  double y_1 = first_agent->y_coord;
  double x_2 = second_agent->x_coord;
  double y_2 = second_agent->y_coord;

  if (fabs(x_1 - x_2) >= (BOARD_LENGTH - 1)) {
    if (x_1 < x_2) {
      x_1 = x_1 + BOARD_LENGTH;
    } else {
      x_2 = x_2 + BOARD_LENGTH;
    }
  }
  if (fabs(y_1 - y_2) >= (BOARD_LENGTH - 1)) {
    if (y_1 < y_2) {
      y_1 = y_1 + BOARD_LENGTH;
    } else {
      y_2 = y_2 + BOARD_LENGTH;
    }
  }
  if (calc_dist(x_1, y_1, x_2, y_2) <= 1) {
    return 1;
  } else {
    return 0;
  }
}

//Calculate distance between two particles
double calc_dist(double x_1, double y_1, double x_2, double y_2) {
  double x_diff = (x_1 - x_2);
  double y_diff = (y_1 - y_2);
  return sqrt((x_diff * x_diff) + (y_diff * y_diff));
}

void update_position(agent *next_agent) {
  next_agent->x_coord = next_agent->x_coord + next_agent->x_dir;
  next_agent->y_coord = next_agent->y_coord + next_agent->y_dir;

  //Account for periodic boundaries
  if (next_agent->x_coord >= BOARD_LENGTH) {
    next_agent->x_coord = next_agent->x_coord - BOARD_LENGTH;
  } else if (next_agent->x_coord < 0) {
    next_agent->x_coord = next_agent->x_coord + BOARD_LENGTH;
  }

  if (next_agent->y_coord >= BOARD_LENGTH) {
    next_agent->y_coord = next_agent->y_coord - BOARD_LENGTH;
  } else if (next_agent->y_coord < 0) {
    next_agent->y_coord = next_agent->y_coord + BOARD_LENGTH;
  }
}

void clear_squares(square **board) {
  for (int i = 0; i < BOARD_LENGTH; i++) {
    for (int j = 0; j < BOARD_LENGTH; j++) {
      board[i][j].index = -1;
    }
  }
}

void update_squares(agent *agents) {
  for (int i = 0; i < AGENT_NUMBER; i++) {
    int x_coord = (int) floor(agents[i].x_coord);
    int y_coord = (int) floor(agents[i].y_coord);
    BOARD[x_coord][y_coord].index = BOARD[x_coord][y_coord].index + 1; //Incrementing count of agents on square
    int index = BOARD[x_coord][y_coord].index;
    BOARD[x_coord][y_coord].active_agents[index] = agents[i].agent_num; //Add agent to list of agents 
  }
}

void print_board() {
  for (int i = 0; i < BOARD_LENGTH; i++) {
    for (int j = 0; j < BOARD_LENGTH; j++) {
      printf("%d ", BOARD[i][j].index + 1);
    }
    printf("\n");
  }
  printf("\n");
}

void print_agents(agent *agents) {
  double order = order_param(agents);
  for (int i = 0; i < AGENT_NUMBER; i++) {
    printf("%f,%f,%f,%f,%f\n", agents[i].x_coord, agents[i].y_coord, agents[i].x_dir, agents[i].y_dir, order); //Print relevant information from previous iteration
  }
  printf("\n");
}

int main(int argc, char **argv) {
  int iterations;

  if (argc < 4) {
    iterations = 1000;
  } else {
    //Store the command line arguments
    BOARD_LENGTH = atoi(argv[1]);
    AGENT_NUMBER = atoi(argv[2]);
    sscanf(argv[3], "%lf", &NOISE);
    iterations = atoi(argv[4]);
  }

  //Initialize the list of agents
  agent *agent_list = malloc(AGENT_NUMBER * sizeof(agent));

  //Seed the random number generator 
  srand(time(NULL));
 
  //Print the relevant system parameters to the beginning of the CSV
  printf("%d\n%d\n%d\n%f\n", AGENT_NUMBER, BOARD_LENGTH, iterations, NOISE);
  
  //Initialize all the agents
  for (int i = 0; i < AGENT_NUMBER; i++) {
    double x_comp = ((double) rand() / (double) RAND_MAX);
    double y_comp = ((double) rand() / (double) RAND_MAX);
    int bool = rand();
    if (bool % 2) {
      x_comp = x_comp * -1;
    } 
    bool = rand();
    if (bool % 2) {
      y_comp = y_comp * -1;
    }
    double magnitude = sqrt((x_comp * x_comp) + (y_comp * y_comp));

    (agent_list+i)->agent_num  = i; 
    (agent_list+i)->x_coord = ((double) rand() / (double) RAND_MAX) * BOARD_LENGTH;
    (agent_list+i)->y_coord = ((double) rand() / (double) RAND_MAX) * BOARD_LENGTH;
    (agent_list+i)->x_dir = (x_comp / magnitude);
    (agent_list+i)->y_dir = (y_comp / magnitude);
  }

  //Initialize all the squares
  BOARD = malloc(BOARD_LENGTH * sizeof(square *));
  for (int i = 0; i < BOARD_LENGTH; i++) {
    BOARD[i] = malloc(BOARD_LENGTH * sizeof(square));
  }

  for (int i = 0; i < BOARD_LENGTH; i++) {
    for (int j = 0; j < BOARD_LENGTH; j++) {
      BOARD[i][j].active_agents = malloc(AGENT_NUMBER * sizeof(int));
      BOARD[i][j].index = -1; //Square indexing system. BOARD[i][j].index is equal to one less than the number agents on that square.
    }
  }

  clear_squares(BOARD);
  update_squares(agent_list); //reallocate agents to their respective squares
  print_agents(agent_list); //print the last iteration to ouptu file

  //Step through iterations of the algorithm
  for (int i = 0; i < iterations; i++) {
    for (int j = 0; j < AGENT_NUMBER; j++) {
      update_orientation(&agent_list[j], agent_list);
      //update_position(&agent_list[j]);
    }
    for (int j = 0; j < AGENT_NUMBER; j++) {
      update_position(&agent_list[j]);
    }
    clear_squares(BOARD);
    update_squares(agent_list);
    print_agents(agent_list);

    COUNTER++;
  }

  //Free agents
  free(agent_list);

  //Free board
  for (int i = 0; i < BOARD_LENGTH; i++) {
    for (int j = 0; j < BOARD_LENGTH; j++) {
      free(BOARD[i][j].active_agents);
    }
    free(BOARD[i]);
  }
  free(BOARD);
  
}
