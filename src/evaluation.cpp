#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>

#include "../include/evaluation.h"
#include "../include/aggregation.h"
#include "../include/solver.h"
#include "../include/alignments.h"
#include "../include/helperfunctions.h"
#include "../include/tec.h"
#include "../include/DFO.h"

using namespace std;

double computeBaselineCost(const vector<Flexoffer> &flexOffers, const vector<double> &spotPrices)
{
    double total_cost = 0.0;
    for (const auto &fo : flexOffers) {
        int duration = fo.get_duration();
        auto profile = fo.get_profile();
        for (int h = 0; h < duration; h++) {
            double avg_power = (profile[h].min_power + profile[h].max_power) / 2.0;
            double price = spotPrices[h + fo.get_est_hour()];
            total_cost += avg_power * price;
        }
    }
    return total_cost;
}

double computeBaselineCost(const vector<Tec_flexoffer> &flexOffers, const vector<double> &spotPrices)
{
    double total_cost = 0.0;
    for (const auto &fo : flexOffers) {
        int duration = fo.get_duration();
        auto profile = fo.get_profile();
        for (int h = 0; h < duration; h++) {
            double avg_power = (profile[h].min_power + profile[h].max_power) / 2.0;
            double price = spotPrices[h + fo.get_est_hour()];
            total_cost += avg_power * price;
        }
    }
    return total_cost;
}

double computeBaselineCost(const vector<DFO> &dfos, const vector<double> &spotPrices) {
    double total_cost = 0.0;

    for (const auto &dfo : dfos) {
        int T = dfo.polygons.size();
        int maxT = min(T, static_cast<int>(spotPrices.size()));

        double accumulated_dependency = 0.0;

        for (int t = 0; t < maxT; ++t) {
            const auto &polygon = dfo.polygons[t];

            vector<Point> matching_points = find_or_interpolate_points(polygon.points, accumulated_dependency);
            double min_current_energy = matching_points[0].y;
            double max_current_energy = matching_points[1].y;
            double avg_current_energy = 0.5 * (min_current_energy + max_current_energy);
            total_cost += avg_current_energy * spotPrices[t];
            accumulated_dependency += avg_current_energy;
        }
    }

    return total_cost;
}

double computeAggregatedCost(vector<Flexoffer> flexOffers, int est_threshold, int lst_threshold, int max_group_size, Alignments align, const vector<double> &spotPrices){
    vector<AggregatedFlexOffer> afos = nToMAggregation(flexOffers, est_threshold, lst_threshold, max_group_size, align, spotPrices, 0);
    Solver::solve(afos, spotPrices);

    double total_cost = 0.0;
    for (auto &afo : afos) {
        const auto &sched = afo.get_scheduled_allocation();
        int duration = afo.get_duration();
        for (int t = 0; t < duration; t++) {
            double power = sched[t];
            double price = spotPrices[t + afo.get_aggregated_earliest_hour()];
            total_cost += power * price;
        }
    }
    return total_cost;
}

double computeAggregatedCost(vector<Tec_flexoffer> flexOffers, int est_threshold, int lst_threshold, int max_group_size, Alignments align, const vector<double> &spotPrices){
    vector<AggregatedFlexOffer> afos = nToMAggregation(flexOffers, est_threshold, lst_threshold, max_group_size, align, spotPrices, 1);
    Solver::solve_tec(afos, spotPrices);

    double total_cost = 0.0;
    for (auto &afo : afos) {
        const auto &sched = afo.get_scheduled_allocation();
        int duration = afo.get_duration();
        for (int t = 0; t < duration; t++) {
            double power = sched[t];
            double price = spotPrices[t + afo.get_aggregated_earliest_hour()];
            total_cost += power * price;
        }
    }
    return total_cost;
}



double computeAggregatedCost(vector<DFO> dfos, const vector<double> &spotPrices, int max_group_size){
    double total_cost = 0.0;
    int usedOffers = 0;


    while (usedOffers < (int)dfos.size())
    {
        int remain = dfos.size() - usedOffers;
        int this_chunk_size = min(remain, max_group_size);

        vector<DFO> dfo_group(dfos.begin() + usedOffers, dfos.begin() + usedOffers + this_chunk_size);

        usedOffers += this_chunk_size;

        double epsilon1 = 1.0;
        double epsilon2 = 1.0;
        DFO AFO = aggnto1(dfo_group, 5, epsilon1, epsilon2);

        //AFO.print_dfo();

        vector<double> schedule = Solver::DFO_Optimization(AFO, spotPrices);
        for (size_t t = 0; t < schedule.size() && t < spotPrices.size(); ++t) {
            total_cost += schedule[t] * spotPrices[t];
        }
    }
    return total_cost;
}


void runAggregationScenarios(const vector<Flexoffer> &normalOffers, const vector<Tec_flexoffer> &tecOffers, const vector<DFO> &dfos, const vector<double> &spotPrices){
    
    auto scenarios = generateScenarioMatrix(); 
    string csvFile = "../data/economic_savings.csv";
    ofstream file(csvFile);
    if (!file.is_open()) {
        cerr << "Error: could not open " << csvFile << endl;
        return;
    }

    file << "scenario_id,aggregator_type,alignment,est_threshold,lst_threshold,max_group_size,baseline_cost,aggregated_cost,savings,scenario_time,NrOfFlexOffers\n";

    int scenario_id = 1;
    for (auto &s : scenarios) {
        double baseline=0.0;
        double agg_cost=0.0;

        int scenarioSize = 0;
        if (s.aggregator_type == 0) {
            scenarioSize = min(s.usedOffers, static_cast<int>(normalOffers.size()));
        } else if (s.aggregator_type == 1) {
            scenarioSize = min(s.usedOffers, static_cast<int>(tecOffers.size()));
        } else if (s.aggregator_type == 2) {
            scenarioSize = min(s.usedOffers, static_cast<int>(dfos.size()));
        }

        if (scenarioSize == 0) {
            cout << "Skipping scenario " << scenario_id << " [type=" << s.aggregator_type
                 << "] due to unavailable flex-offers." << endl;
            continue;
        }

        auto t_start = chrono::steady_clock::now();
        if (s.aggregator_type == 0) { //thsi is for normal FOs
            vector<Flexoffer> subNormal(normalOffers.begin(), normalOffers.begin() + scenarioSize);
            baseline = computeBaselineCost(subNormal, spotPrices);
            agg_cost = computeAggregatedCost(subNormal, s.est_threshold, s.lst_threshold, s.max_group_size, s.align, spotPrices);
        }
        else if (s.aggregator_type == 1) { // this is for tec FOs
            vector<Tec_flexoffer> subTec(tecOffers.begin(), tecOffers.begin() + scenarioSize);
            baseline = computeBaselineCost(subTec, spotPrices);
            agg_cost = computeAggregatedCost(subTec, s.est_threshold, s.lst_threshold, s.max_group_size, s.align, spotPrices);
        }
        else if (s.aggregator_type == 2){ // DFO
            vector<DFO> subDFOs(dfos.begin(), dfos.begin() + scenarioSize);
            baseline = computeBaselineCost(subDFOs, spotPrices);
            agg_cost = computeAggregatedCost(subDFOs, spotPrices, s.max_group_size);
        }

        auto t_end = chrono::steady_clock::now();
        double savings = baseline - agg_cost;

        double scenarioTimeSec = chrono::duration<double>(t_end - t_start).count();

        file << scenario_id << "," << s.aggregator_type << "," << static_cast<int>(s.align) << "," << s.est_threshold << "," << s.lst_threshold << ","
        << s.max_group_size << "," << baseline << "," << agg_cost << "," << savings << "," << scenarioTimeSec << "," << scenarioSize << "\n";

        cout << "Scenario " << scenario_id << " [type=" << s.aggregator_type << " alignment= " << static_cast<int>(s.align) << "] => " << "est=" << s.est_threshold << ", lst=" << s.lst_threshold
        << ", maxG=" << s.max_group_size << "\n" << "   baseline=" << baseline  << ", aggregated_cost=" << agg_cost  << ", savings=" << savings  << ", NrOfFos " <<  scenarioSize << ", scenario_time=" << scenarioTimeSec << "s\n\n";
        scenario_id++;
    }

    file.close();
}

vector<AggScenario> generateScenarioMatrix() {

    vector<AggScenario> scenarios;
    vector<int> aggrTypes = {0, 1, 2};
    vector<Alignments> aligns = {
        Alignments::start,
        Alignments::balance,
        Alignments::price,
    };

    vector<int> thresholds = {2, 4, 6}; 
    vector<int> groupSizes = {5, 10, 50};
    vector<int> nOffersVec = {10, 50, 100, 200};

    for (int at : aggrTypes) {
        for (auto al : aligns) {
            for (int th : thresholds) {
                for (int g : groupSizes) {
                    for (int usedN : nOffersVec) {
                        AggScenario s;
                        s.aggregator_type = at;
                        s.est_threshold = th;
                        s.lst_threshold = th;
                        s.max_group_size = g;
                        s.align = al;
                        s.usedOffers = usedN;
                        scenarios.push_back(s);
                    }
                }
            }
        }
    }
    return scenarios;
}
