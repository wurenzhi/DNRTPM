#include <iostream>
#include <array>
#include <vector>
#include <cassert>
#include <memory>
#include <math.h>
#include <limits>
#include <deque>
#include <random>
#include <queue>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <boost/circular_buffer.hpp>
#include <boost/intrusive/set.hpp>

using namespace std;
constexpr auto inf = numeric_limits<double>::infinity();
constexpr auto inf_int = numeric_limits<int>::max();

using Match = pair<pair<double,int>, pair<int,int>>;
bool Match_compare(Match a, Match b){
    return a.first.first < b.first.first;
}

inline double dist (double x, double y){
    return (x-y)*(x-y);
};

pair<double,double> mean_std(std::vector<double>& samples)
{
    int size = samples.size();

    double sums = 0;
    double sums_square = 0;
    for (int i = 0; i < size; i++)
    {
        sums += samples[i];
        sums_square += samples[i]*samples[i];
    }
    double mean = sums/size;
    return make_pair(mean,sqrt(sums_square / size-mean*mean));
}

pair<vector <double>,vector <double>> prefix_norm(vector <double>& query_ori){
    vector<double > query ;
    bool append_flag = false;
    query.push_back(query_ori[0]);
    for (int j = 1; j < query_ori.size(); ++j) {
        if(query_ori[j]!=query_ori[0]){
            append_flag = true;
        }
        if(append_flag){
            query.push_back(query_ori[j]);
        }
    }
    int n = query.size();
    auto m_sd_pair = mean_std(query);
    double std_all = m_sd_pair.second;
    vector <double> M_sums(1);
    vector <double> SD_sums(1);
    vector <double> n_query(n);
    vector <double> factor_sigma(n);

    for(int i=0;i<n;i++){
        M_sums[M_sums.size()-1]+= query[i];
        SD_sums[SD_sums.size()-1]+=query[i]*query[i];
        M_sums.push_back(M_sums[M_sums.size()-1]);
        SD_sums.push_back(SD_sums[SD_sums.size()-1]);

        double mu = M_sums[i]/float(i+1);
        double sd = sqrt(fabs(SD_sums[i]/float(i+1) - mu*mu));
        if(sd < 1e-10){
            n_query[i]=query[i]-mu;
        } else{
            n_query[i]=(query[i]-mu)/sd;
        }
        factor_sigma[i]=sd/std_all;
    }
    return make_pair(n_query,factor_sigma);
}


vector<Match> DNRTPM(vector <double> stream, vector <double> query, bool best_query,int n_query, bool disjoint, double threshold){
    double epsilon = threshold;
    if(best_query)
        epsilon = inf;
    auto pr = prefix_norm(query);
    auto template_norm = pr.first;
    auto factor_sigma_squre = pr.second;
    for(auto & sigma:factor_sigma_squre){
        sigma = sigma*sigma;
    }

    int n = template_norm.size();
    vector<double> D_recent(n);
    vector<double> D_now(n);
    vector<int> S_recent(n),S_now(n);
    double d_rep = epsilon;
    int J_s = inf_int;
    int J_e = inf_int;
    for(auto & d:D_recent){
        d = inf;
    }
    int check = 0;
    boost::circular_buffer<double> M_sums(10000),SD_sums(10000);
    int sums_start_index = 0;
    priority_queue<Match ,vector<Match >,function<bool (Match,Match)>> best_matches(Match_compare);
//    Match best_match;
    vector<Match> matches;

    for(int j = 0;j < stream.size();j++){
        double x = stream[j];
        if(j == 0) {
            M_sums.push_back(x);
            SD_sums.push_back(x * x);
        } else{
            if(M_sums.full()){
                M_sums.set_capacity(M_sums.size()+10000);
            }
            if(SD_sums.full()){
                SD_sums.set_capacity(SD_sums.size()+10000);
            }
            M_sums.push_back(M_sums[M_sums.size()- 1] + x);
            SD_sums.push_back(SD_sums[SD_sums.size()- 1] + x * x);
        }


        //---------------------------- filling in current column of stwm
        int start_min = inf_int;

        double d_now_i_1,d_recent_i_1,d_recent_i,min_d, mean,std;
        double x_norm = 0.0;
        int t_start;
        for(int i = 0;i<n;i++){

            if(i==0){
                S_now[i] = j;
                D_now[i] = 0;
            } else{

                    if(D_now[i-1]< epsilon){

                        t_start = S_now[i-1];
                        if(j == t_start){
                            x_norm =  0.0;
                        } else{

                            if(t_start == 0){
                                mean = M_sums[j-sums_start_index]/(j+1);
                                std = sqrt(fabs(SD_sums[j-sums_start_index] / (j + 1) - mean*mean));
                            } else{
                                mean = (M_sums[j - sums_start_index] - M_sums[t_start - 1 - sums_start_index]) / (j - t_start + 1);
                                std = sqrt(fabs((SD_sums[j - sums_start_index] - SD_sums[t_start - 1 - sums_start_index]) / (
                                        j - t_start + 1) - mean * mean));
                            }
                            if(std <= 1e-10){
                                x_norm =  (x - mean);
                            } else{
                                x_norm =  (x - mean)/std;
                            }
                        }
                        //x_norm = prefix_norm_t(x,j,S_now[i-1],sums_start_index,M_sums,SD_sums);
                        d_now_i_1 = D_now[i-1] + (x_norm-template_norm[i])*(x_norm-template_norm[i]) * factor_sigma_squre[i];
                    }else{
                        d_now_i_1 = inf;
                    }
                    if( D_recent[i-1]<epsilon){

                        t_start = S_recent[i-1];
                        if(j == t_start){
                            x_norm =  0.0;
                        } else{

                            if(t_start == 0){
                                mean = M_sums[j-sums_start_index]/(j+1);
                                std = sqrt(fabs(SD_sums[j-sums_start_index] / (j + 1) - mean*mean));
                            } else{
                                mean = (M_sums[j - sums_start_index] - M_sums[t_start - 1 - sums_start_index]) / (j - t_start + 1);
                                std = sqrt(fabs((SD_sums[j - sums_start_index] - SD_sums[t_start - 1 - sums_start_index]) / (
                                        j - t_start + 1) - mean * mean));
                            }
                            if(std <= 1e-10){
                                x_norm =  (x - mean);
                            } else{
                                x_norm =  (x - mean)/std;
                            }
                        }

                        //x_norm = prefix_norm_t(x,j,S_recent[i-1],sums_start_index,M_sums,SD_sums);
                        d_recent_i_1 = D_recent[i-1] + (x_norm-template_norm[i])*(x_norm-template_norm[i]) * factor_sigma_squre[i];
                    }else{
                        d_recent_i_1 = inf;
                    }
                    if(D_recent[i]<epsilon){

                        t_start = S_recent[i];
                        if(j == t_start){
                            x_norm =  0.0;
                        } else{

                            if(t_start == 0){
                                mean = M_sums[j-sums_start_index]/(j+1);
                                std = sqrt(fabs(SD_sums[j-sums_start_index] / (j + 1) - mean*mean));
                            } else{
                                mean = (M_sums[j - sums_start_index] - M_sums[t_start - 1 - sums_start_index]) / (j - t_start + 1);
                                std = sqrt(fabs((SD_sums[j - sums_start_index] - SD_sums[t_start - 1 - sums_start_index]) / (
                                        j - t_start + 1) - mean * mean));
                            }
                            if(std <= 1e-10){
                                x_norm =  (x - mean);
                            } else{
                                x_norm =  (x - mean)/std;
                            }
                        }

                        //x_norm = prefix_norm_t(x,j,S_recent[i],sums_start_index,M_sums,SD_sums);
                        d_recent_i = D_recent[i] +  (x_norm-template_norm[i])*(x_norm-template_norm[i]) * factor_sigma_squre[i];
                    }else{
                        d_recent_i = inf;
                    }
                    min_d = min(min(d_now_i_1 , d_recent_i_1), d_recent_i);
                    if(min_d<epsilon){
                        if(min_d == d_now_i_1)
                            S_now[i] = S_now[i - 1];
                        else if(min_d == d_recent_i_1)
                            S_now[i] = S_recent[i - 1];
                        else if(min_d == d_recent_i)
                            S_now[i] = S_recent[i];
                        D_now[i] = min_d;
                    } else{
                        D_now[i] = inf;
                        S_now[i] = j;
                    }
            }
            start_min = min(start_min, S_now[i]);
        }
        while(start_min-1 > sums_start_index){
            M_sums.pop_front();
            SD_sums.pop_front();
            sums_start_index = sums_start_index+1;
        }
        //-------------------------------------------



        if(D_now[n - 1] <= epsilon and  D_now[n - 1] < d_rep){
            d_rep = D_now[n - 1];
            J_s = S_now[n - 1];
            J_e = j + 1;
            if(!disjoint){
                if(matches.empty() or J_s > matches[matches.size()-1].second.second or d_rep < matches[matches.size()-1].first.first){
                    matches.emplace_back(make_pair(make_pair(d_rep,0), make_pair(J_s, J_e)));
                }
            }
        }

        // disjoint query
        for(int i=0;i<n;i++){
            if(D_now[i] >= d_rep or S_now[i] > J_e)
                check = check + 1;
        }
        if(check==n){
            if(best_query){
                if(best_matches.size() < n_query){
                    best_matches.emplace(make_pair(make_pair(d_rep,j-J_e+1), make_pair(J_s, J_e)));
                }
                else if(best_matches.top().first.first > d_rep){
                    best_matches.pop();
                    best_matches.emplace(make_pair(make_pair(d_rep,j-J_e+1), make_pair(J_s, J_e)));
                    epsilon = best_matches.top().first.first;
                }
            }else if(disjoint){
                matches.emplace_back(make_pair(make_pair(d_rep,j-J_e+1), make_pair(J_s, J_e)));
            }

            for(int i=0;i<n;i++) {
                if(D_now[i] >= d_rep and S_now[i] < J_e)
                    D_now[i] = inf;
            }

            d_rep = inf;
            J_s = inf_int;
            J_e = inf_int;
        }
        check = 0;
        swap(D_recent,D_now);
        swap(S_recent,S_now);
    }
    if(best_query){
        matches.clear();
        while (!best_matches.empty()){
            matches.emplace_back(best_matches.top());
            best_matches.pop();
        }
    }
    return matches;
}
