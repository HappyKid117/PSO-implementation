#include<bits/stdc++.h>
#include<random>
using namespace std;


double upperBound = 5.0;
double lowerBound = -5.0;
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> pos_distr(lowerBound, upperBound);
uniform_real_distribution<> coeff_distr(0, 1);


//----------------------------------------------------
// Rastrigin Function
// double fitness(vector<double> position){
//     int dim = position.size();
//     double fitnessValue = 10*dim;
//     for(int i=0; i<dim; i++){
//         double x = position[i];
//         fitnessValue += (x*x) - (10*cos(2*M_PI*x));
//     }
//     return fitnessValue;
// }

// Finding Food at (x,y)
// double fitness(vector<double> position){
//     vector<double> foodPosition = {700.0, 700.0};
//     int dim = position.size();
//     double fitnessValue = 0;
//     for(int i=0; i<dim; i++){
//         fitnessValue += abs(position[i] - foodPosition[i]);
//     }
//     return fitnessValue;
// }

// Rosenbrock function
// double fitness(vector<double> position){
//     int dim = position.size();
//     double fitnessValue = 0;
//     for(int i=0; i<dim-1; i++){
//         double xi = position[i];
//         double xi1 = position[i+1];
//         fitnessValue += 100*pow((xi1 - (xi*xi)),2) + pow(1-xi, 2);
//     }
//     return fitnessValue;
// }


// Styblinski–Tang function
double fitness(vector<double> position){
    int dim = position.size();
    double fitnessValue = 10*dim;
    for(int i=0; i<dim; i++){
        double x = position[i];
        fitnessValue += pow(x,4) - 16*pow(x,2) + 5*x;
    }
    return fitnessValue/2;
}

// sin(x)
// double fitness(vector<double> position){
//     int dim = position.size();
//     double fitnessValue = 0;
//     for(int i=0; i<dim-1; i++){
//         fitnessValue += sin(position[i]);
//     }
//     return fitnessValue;
// }

//-----------------------------------------------------

class Particle{
    public:
        vector<double> position;
        double currentFitness;
        vector<double> pbest;
        double bestFitness;
        vector<double> velocity;


        Particle(int dim){
            for(int i=0; i<dim; i++){
                position.push_back(pos_distr(gen));
                velocity.push_back(pos_distr(gen));
            }
            currentFitness = fitness(position);
            pbest = position;
            bestFitness = currentFitness;
        }

        void printPosition(){
            for(int i=0; i<position.size(); i++){
                cout<<position[i]<<" ";
            }
            cout<<endl;
        }

};

vector<vector<Particle>> history;

vector<double> pso(int epochs, int n, int dim){
    double w = 0.729;
    double c1 = 1.49445;
    double c2 = 1.49445;
    
    vector<Particle> Swarm;
    vector<double> gbest(dim);
    double bestFitness = INT_MAX;

    for(int i=0; i<n; i++){
        Swarm.push_back(Particle(dim));
    }

    for(int i=0; i<n; i++){
        if(Swarm[i].currentFitness < bestFitness){
            bestFitness = Swarm[i].currentFitness;
            gbest = Swarm[i].position;
        }
    }

    for(int t=0; t<epochs; t++){
        history.push_back(Swarm);
        for(int i=0; i<n; i++){
            double r1 = coeff_distr(gen);
            double r2 = coeff_distr(gen);

            for(int j=0; j<dim; j++){
                Swarm[i].velocity[j] = (w*Swarm[i].velocity[j])
                                     + (c1*r1*(Swarm[i].pbest[j] - Swarm[i].position[j]))
                                     + (c2*r2*(gbest[j] - Swarm[i].position[j]));
                
                Swarm[i].velocity[j] = min(upperBound, Swarm[i].velocity[j]);
                Swarm[i].velocity[j] = max(lowerBound, Swarm[i].velocity[j]);

                Swarm[i].position[j] += Swarm[i].velocity[j];
            }

            Swarm[i].currentFitness = fitness(Swarm[i].position);

            if(Swarm[i].currentFitness < Swarm[i].bestFitness){
                Swarm[i].bestFitness = Swarm[i].currentFitness;
                Swarm[i].pbest = Swarm[i].position; 
            }

            if(Swarm[i].currentFitness < bestFitness){
                bestFitness = Swarm[i].currentFitness;
                gbest = Swarm[i].position;
            }
        }
    }

    return gbest;
}

int main(){
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    vector<double> answer;

    int number_of_iterations = 200;
    int number_of_particles = 100;
    int dimensions = 5;

    history.clear();
    answer = pso(number_of_iterations, number_of_particles, dimensions);

    cout<<"Best Position : ";
    for(int i=0; i<dimensions; i++){
        cout<<answer[i]<<" ";
    }
    cout<<endl;
    cout<<"Fitness : "<<fitness(answer)<<endl;

    ofstream file;
    file.open("history.txt");
    file<<number_of_iterations<<" "<<number_of_particles<<" "<<dimensions<<" "<<lowerBound<<" "<<upperBound<<endl;
    for(int t=0; t<number_of_iterations; t++){
        for(int i=0; i<number_of_particles; i++){
            for(int j=0; j<dimensions; j++){
                file<<history[t][i].position[j]<<" ";
            }
            file<<endl;
        }
    }
    cout<<"History printed in \"history.txt\""<<endl;

}