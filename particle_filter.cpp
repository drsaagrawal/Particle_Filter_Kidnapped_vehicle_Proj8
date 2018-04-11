/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <limits>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    num_particles = 15;
    default_random_engine gen;
    // This line creates a normal (Gaussian) distribution for x.
	normal_distribution<double> dist_x(x, std[0]);

	// TODO: Create normal distributions for y and theta.
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

    for (int i=0;i<num_particles;i++)
    {
        double sample_x, sample_y, sample_theta;
        Particle p;
        // where "gen" is the random engine initialized earlier.
		sample_x = dist_x(gen);
        sample_y = dist_y(gen);
		sample_theta = dist_theta(gen);
		p.id = i;
		p.x = sample_x;
		p.y = sample_y;
		p.theta = sample_theta;
		p.weight = 1;
		particles.push_back(p);
		weights.push_back(1);
    }
    is_initialized = true;
    //cout<<"Initialization out"<<endl;
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    //cout<<"Prediction in"<<endl;
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;
	for (int i=0;i<num_particles;i++)
    {
        double sample_x, sample_y, sample_theta,x,y,theta;
       if (yaw_rate == 0)
        {
             x = particles[i].x+ (velocity*delta_t*cos(particles[i].theta));
             y = particles[i].y+ (velocity*delta_t*sin(particles[i].theta));
             theta = particles[i].theta;
        }
        else{
             x = particles[i].x+((velocity/yaw_rate)*(sin(particles[i].theta+(yaw_rate*delta_t))-sin(particles[i].theta)));
             y = particles[i].y+((velocity/yaw_rate)*(-cos(particles[i].theta+(yaw_rate*delta_t))+cos(particles[i].theta)));
             theta = particles[i].theta+(yaw_rate*delta_t);
        }

        // This line creates a normal (Gaussian) distribution for x.
        normal_distribution<double> dist_x(x, std_pos[0]);
        // TODO: Create normal distributions for y and theta.
        normal_distribution<double> dist_y(y, std_pos[1]);
        normal_distribution<double> dist_theta(theta, std_pos[2]);
        sample_x = dist_x(gen);
        sample_y = dist_y(gen);
		sample_theta = dist_theta(gen);
		particles[i].x = (sample_x);
        particles[i].y = (sample_y);
        particles[i].theta = (sample_theta);
    }
    //cout<<"Prediction out"<<endl;

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
    //cout<<"Data Association in"<<endl;
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

	//predicted -> map     observation -> map transform
	for (int i=0;i<observations.size();i++)
    {
        double obs_x = observations[i].x;    // Obs in map coordinate
        double obs_y = observations[i].y;
        double dist = std::numeric_limits<double>::infinity();

        for (int j=0;j<predicted.size();j++)
        {
            double map_x = predicted[j].x;     // map landmarks x
            double map_y = predicted[j].y;
            double temp_dist = sqrt((obs_x-map_x)*(obs_x-map_x) + ((obs_y-map_y)*(obs_y-map_y)));
            if (temp_dist<dist)
            {
                dist = temp_dist;
                observations[i].id = predicted[j].id;
            }

        }
    }
    //cout<<"Data Association out"<<endl;

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
    //cout<<"Update weight in"<<endl;

    std::vector<LandmarkObs> map_land;
    LandmarkObs l_map;
    for (int i=0; i < map_landmarks.landmark_list.size();i++)
    {
        l_map.id = map_landmarks.landmark_list[i].id_i;
        l_map.x = map_landmarks.landmark_list[i].x_f;
        l_map.y = map_landmarks.landmark_list[i].y_f;
        map_land.push_back(l_map);
    }


	for (int i =0; i < num_particles;i++)
    {
        std::vector<LandmarkObs> temp_observations;
        LandmarkObs temp;
        //temp_observations = observations;
        for (int j=0;j<observations.size();j++)
        {
            //temp.id = observations[j].id;
            LandmarkObs temp;
            temp.x = particles[i].x + (cos(particles[i].theta)*observations[j].x)-(sin(particles[i].theta)*observations[j].y);
            temp.y = particles[i].y + (cos(particles[i].theta)*observations[j].y)+(sin(particles[i].theta)*observations[j].x);
            temp_observations.push_back(temp);
        }

        dataAssociation(map_land,temp_observations);
        std::vector<int> temp_associations;
        std::vector<double> temp_sense_x;
        std::vector<double> temp_sense_y;
        //cout<<"Hi"<<temp_observations[0].id;
        for (int j=0;j<temp_observations.size();j++)
        {
            temp_associations.push_back(temp_observations[j].id);
            temp_sense_x.push_back(temp_observations[j].x);
            temp_sense_y.push_back(temp_observations[j].y);
        }
        SetAssociations(particles[i],temp_associations,temp_sense_x,temp_sense_y);
        //cout<<"Back from setassociation"<<endl;
        double sig_x= std_landmark[0];
        //cout<<"1"<<endl;
        double sig_y= std_landmark[1];
        //cout<<"2"<<endl;
        double temp_weight = 1.0;
        //cout<<"3"<<endl;

        //calculationg weight

        for (int j=0;j<temp_observations.size();j++)
        {
            //cout << "Calculate for1 in"<<endl;

            for (int k=0;k<map_land.size();k++)
            {
                //cout << "Calculate for2 in"<<endl;
                if(temp_observations[j].id == map_land[k].id)
                {
                    double x_obs = temp_observations[j].x;
                    double y_obs= temp_observations[j].y;
                    double mu_x= map_land[k].x;
                    double mu_y= map_land[k].y;
                    // calculate normalization term
                    double gauss_norm= (1/(2 * M_PI * sig_x * sig_y));
                    // calculate exponent
                    double exponent= (pow((x_obs - mu_x),2))/(2 * pow(sig_x,2)) + (pow((y_obs - mu_y),2))/(2 * pow(sig_y,2));
                    //calculate weight using normalization terms and exponent
                    double weight= gauss_norm * exp(-exponent);
                    temp_weight *=weight;
                    break ;
                }
            }
        }

        // Updating weighjt
        particles[i].weight = temp_weight;
        weights[i]=temp_weight;

    }
    //cout<<"Update weight out"<<endl;
}

void ParticleFilter::resample() {
    //cout<<"Resample in"<<endl;
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::default_random_engine generator;
	std::discrete_distribution<int> distribution {weights.begin(),weights.end()};
	vector<Particle> p2;
    double new_weight_sum = 0;

    for(int i=0;i<num_particles;i++)
    {
       p2.push_back(particles[distribution(generator)]);
       weights[i]  = p2[i].weight;
       new_weight_sum += weights[i];
    }
   //normalize
    for(int i=0;i<num_particles;i++)
    {
       weights[i] /= new_weight_sum;
       p2[i].weight = weights[i];
    }
    particles = p2;
    //cout<<"Resample out"<<endl;
}

void ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //cout<<"Setassociation in"<<endl;
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates
    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
