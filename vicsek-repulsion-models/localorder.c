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
  double x_dir_new;
  double y_dir;
  double y_dir_new;
};

struct square {
  int *active_agents;
  int index;
};

double in_radius(agent *first_agent, agent *second_agent);
double calc_dist(double x_1, double y_1, double x_2, double y_2);
 
int BOARD_LENGTH = 140;
int AGENT_NUMBER = 10000;
double ORDER_PARAM = 0.0;
double LOCAL_ORDER = 0.0;
double DISPLACEMENT = 0.0;
double ALIGNMENT = 0.5;
double REPULSION = 0.5;
double STRESS = 0.0;
int COUNTER = 0;
square **BOARD;

//Calculate the order parameter based on the orientation of the agents
double order_param(agent *agent_list) {
  double avg_x = 0.0;
  double avg_y = 0.0;
  for (int i = 0; i < AGENT_NUMBER; i++) {
    double magnitude = sqrt((agent_list[i].x_dir * agent_list[i].x_dir) + (agent_list[i].y_dir * agent_list[i].y_dir));
    avg_x = avg_x + (agent_list[i].x_dir / magnitude);
    avg_y = avg_y + (agent_list[i].y_dir / magnitude);
  }
  double total = sqrt(pow(avg_x/(double) AGENT_NUMBER, 2) + pow(avg_y / (double) AGENT_NUMBER, 2));
  return total;
}

//Update the orientaton of the agents to align with its neighbors
double update_orientation_2(agent *next_agent, agent *agent_list) {
  int x_cell = (int) floor(next_agent->x_coord);
  int y_cell = (int) floor(next_agent->y_coord);
  int neighbor_count = 0;
  double align_x = 0.0;
  double align_y = 0.0;
  double repulse_x = 0.0;
  double repulse_y = 0.0;
  double stress = 0.0;

  for (int x = -3; x <= 3; x++) {
    for (int y = -3; y <= 3; y++) {
      int x_tmp = x_cell + x;
      int y_tmp = y_cell + y;
      
      if (x_tmp < 0 || x_tmp >= BOARD_LENGTH) {
	       continue;
      } 
      if (y_tmp < 0) {
        y_tmp = y_tmp + BOARD_LENGTH;
      }
      if (y_tmp >= BOARD_LENGTH) {
        y_tmp = y_tmp - BOARD_LENGTH;
      }
      for (int i = 0; i <= BOARD[x_tmp][y_tmp].index; i++) {
	       if (!(BOARD[x_tmp][y_tmp].active_agents[i] == next_agent->agent_num)) {
	         int neighbor_number = BOARD[x_tmp][y_tmp].active_agents[i];
           double distance = in_radius(next_agent, &agent_list[neighbor_number]);
           double x_comp = 0.0;
           double y_comp = 0.0;

	        if (distance <= 2.5) {
	          neighbor_count++;
            double neigh_magnitude = sqrt(pow(agent_list[neighbor_number].x_dir, 2) + pow(agent_list[neighbor_number].y_dir, 2));
	          align_x = align_x + (agent_list[neighbor_number].x_dir / neigh_magnitude);
	          align_y = align_y + (agent_list[neighbor_number].y_dir / neigh_magnitude);
          }
          if (distance <= 2) {
            x_comp = next_agent->x_coord - agent_list[neighbor_number].x_coord;
            if (fabs(agent_list[neighbor_number].y_coord - next_agent->y_coord) >= BOARD_LENGTH - 2) {
              y_comp = BOARD_LENGTH - (fabs(agent_list[neighbor_number].y_coord - next_agent->y_coord));
              if (next_agent->y_coord > agent_list[neighbor_number].y_coord) {
                y_comp = y_comp * -1;
              }
            } else {
                y_comp = next_agent->y_coord - agent_list[neighbor_number].y_coord;
            }
            
            double magnitude = sqrt((x_comp * x_comp) + (y_comp * y_comp));
            repulse_x = repulse_x + (x_comp / magnitude);
            repulse_y = repulse_y + (y_comp / magnitude);
            stress = (2 - distance) / 2;
          }

	       }
	     } 
      }

    }

  if (neighbor_count > 0) {
    double magnitude_dir = sqrt((next_agent->x_dir * next_agent->x_dir) + (next_agent->y_dir * next_agent->y_dir));
    double magnitude_align = sqrt((align_x * align_x) + (align_y * align_y));
    double projection = fabs(((next_agent->x_dir / magnitude_dir) * (align_x / magnitude_align)) + ((next_agent->y_dir / magnitude_dir) * (align_y / magnitude_align)));
    LOCAL_ORDER += projection;

    double magnitude_repulse = 1.0;
    next_agent->x_dir_new = ALIGNMENT * (align_x / magnitude_align) + REPULSION * (repulse_x / magnitude_repulse);
    next_agent->y_dir_new = ALIGNMENT * (align_y / magnitude_align) + REPULSION * (repulse_y / magnitude_repulse);
  } else {
    next_agent->x_dir_new = next_agent->x_dir;
    next_agent->y_dir_new = next_agent->y_dir;
  }
  return stress;
}

double in_radius(agent *first_agent, agent *second_agent) {
  double x_1 = first_agent->x_coord;
  double y_1 = first_agent->y_coord;
  double x_2 = second_agent->x_coord;
  double y_2 = second_agent->y_coord;

  if (fabs(y_1 - y_2) >= (BOARD_LENGTH - 3)) {
    if (y_1 < y_2) {
      y_1 = y_1 + BOARD_LENGTH;
    } else {
      y_2 = y_2 + BOARD_LENGTH;
    }
  } 
  return calc_dist(x_1, y_1, x_2, y_2);
}

double calc_dist(double x_1, double y_1, double x_2, double y_2) {
  double x_diff = (x_1 - x_2);
  double y_diff = (y_1 - y_2);
  return sqrt((x_diff * x_diff) + (y_diff * y_diff));
}

void update_orientation(agent *next_agent, agent *agent_list) {
  double new_x = 0.0;
  double new_y = 0.0;
  for (int i = 0; i < AGENT_NUMBER; i++) {
    if (!(i == next_agent->agent_num)) {
      new_x = new_x + agent_list[i].x_dir;
      new_y = new_y + agent_list[i].y_dir;
    }
  }
  
  double current_angle = atan2(new_y, new_x);
  next_agent->x_dir = cos(current_angle);
  next_agent->y_dir = sin(current_angle);
}

void update_position(agent *next_agent) {

  if (next_agent->x_coord + next_agent->x_dir >= BOARD_LENGTH || next_agent->x_coord + next_agent->x_dir < 0) {
    next_agent->x_dir = next_agent->x_dir * -1.0;
  } 
  /*if (next_agent->x_coord < 0.5 || next_agent->x_coord > BOARD_LENGTH - 0.5) {
    next_agent->x_dir_new = next_agent->x_dir * -1.0;
  }*/
  next_agent->x_coord = next_agent->x_coord + next_agent->x_dir;
  next_agent->y_coord = next_agent->y_coord + next_agent->y_dir;

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
    BOARD[x_coord][y_coord].index = BOARD[x_coord][y_coord].index + 1;
    int index = BOARD[x_coord][y_coord].index;
    BOARD[x_coord][y_coord].active_agents[index] = agents[i].agent_num;
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
    printf("%f,%f,%f,%f,%f\n", agents[i].x_coord, agents[i].y_coord, order, STRESS, LOCAL_ORDER / AGENT_NUMBER);
  }
  printf("\n");
}

int main(int argc, char **argv) {
  if (argc < 5) {
    printf("5 arguments expected: (1) Board size (2) Agent number (3) Alignment weight (4) Repulsion weight (5) Iterations\n");
    exit(0);
  }
  //Store the command line arguments
  BOARD_LENGTH = atoi(argv[1]);
  AGENT_NUMBER = atoi(argv[2]);
  sscanf(argv[3], "%lf", &ALIGNMENT);
  sscanf(argv[4], "%lf", &REPULSION);
  int iterations = atoi(argv[5]);
  agent *agent_list = malloc(AGENT_NUMBER * sizeof(agent));

  //Seed the random number generator 
  srand(time(NULL));
 
  printf("%d\n%d\n%d\n", AGENT_NUMBER, BOARD_LENGTH, iterations);
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
    (agent_list+i)->x_dir_new = 0.0;
    (agent_list+i)->y_dir = (y_comp / magnitude);
    (agent_list+i)->y_dir_new = 0.0;
  }

  //Initialize all the squares
  BOARD = malloc(BOARD_LENGTH * sizeof(square *));
  for (int i = 0; i < BOARD_LENGTH; i++) {
    BOARD[i] = malloc(BOARD_LENGTH * sizeof(square));
  }

  for (int i = 0; i < BOARD_LENGTH; i++) {
    for (int j = 0; j < BOARD_LENGTH; j++) {
      BOARD[i][j].active_agents = malloc(AGENT_NUMBER * sizeof(int));
      BOARD[i][j].index = -1;
    }
  }

  clear_squares(BOARD);
  update_squares(agent_list);
  print_agents(agent_list);

  //double STRESS_TOT = 0.0;

  //Step through iterations of the algorithm
  for (int i = 0; i < iterations; i++) {
    if (i == iterations / 2) {
      //REPULSION = 0.025;
    }
    STRESS = 0.0;
    for (int j = 0; j < AGENT_NUMBER; j++) {
      STRESS += update_orientation_2(&agent_list[j], agent_list);
    }
    for (int j = 0; j < AGENT_NUMBER; j++) {
      agent_list[j].x_dir = agent_list[j].x_dir_new;
      agent_list[j].y_dir = agent_list[j].y_dir_new;
      //printf("My x-coordinate is: %f, my x-dir is: %f\nm", agent_list[j].x_coord, agent_list[j].x_dir);
    }
    for (int j = 0; j < AGENT_NUMBER; j++) {
      update_position(&agent_list[j]);
    }


    STRESS = STRESS / (2.0 * AGENT_NUMBER);
    //STRESS_TOT += STRESS;
    //STRESS = STRESS_TOT / (double) ++COUNTER;

    clear_squares(BOARD);
    update_squares(agent_list);
    print_agents(agent_list);

    LOCAL_ORDER = 0.0;
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
