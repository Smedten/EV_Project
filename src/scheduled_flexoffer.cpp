#include <ctime>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "../include/scheduled_flexoffer.h"
#include "../include/flexoffer.h"
#include "../include/aggregation.h"
using namespace std;

//Constructor
ScheduledFlexOffer::ScheduledFlexOffer(AggregatedFlexOffer& AFO)
    : aggregated_offer(AFO){

    const auto &aggregated_profile = aggregated_offer.get_aggregated_profile();
    scheduled_profile.resize(aggregated_profile.size());
    for (size_t i = 0; i < aggregated_profile.size(); i++) {
        scheduled_profile[i] = (aggregated_profile[i].min_power +
                                aggregated_profile[i].max_power) / 2.0;
    }
}

//Getters
AggregatedFlexOffer& ScheduledFlexOffer::get_aggregated_offer() {return aggregated_offer;};
int ScheduledFlexOffer::get_offer_id() {return offer_id;};
vector<double> ScheduledFlexOffer::get_scheduled_profile() {return scheduled_profile;};

//Setters
void ScheduledFlexOffer::set_aggregated_offer(AggregatedFlexOffer& value) {aggregated_offer = value;};
void ScheduledFlexOffer::set_offer_id(int value) {offer_id = value;};
void ScheduledFlexOffer::set_scheduled_profile(vector<double> value) {scheduled_profile = value;};


void ScheduledFlexOffer::print_schedule() {
    std::cout << "Scheduled allocation: " << std::endl;
    for (size_t i = 0; i < aggregated_offer.get_aggregated_profile().size(); i++) {
        if (scheduled_profile[i] > 0){
            std::cout << "  Hour " << i << ": Scheduled Power = " 
                    << fixed << setprecision(2) << scheduled_profile[i] << " kW" << endl;
        }
    }
};


// Implementing 1-To-N Disaggregation without the profile alignment vector
void ScheduledFlexOffer::n_to_1_disaggregation(vector<Flexoffer> &F, AggregatedFlexOffer &fa) {
    // Step 1: Loop through each profile element of the aggregated flex-offer (fa)
    const auto aggregated_profile = fa.get_aggregated_profile();
    const size_t profile_size = aggregated_profile.size();

    vector<double> s_min(profile_size, 0.0);
    vector<double> s_max(profile_size, 0.0);
    vector<double> sx(profile_size, 0.0);
    vector<double> fraction(profile_size, 0.0);

    for (size_t i = 0; i < profile_size; i++) {
        s_min[i] = aggregated_profile[i].min_power;
        s_max[i] = aggregated_profile[i].max_power;
        if (i < scheduled_profile.size()) {
            sx[i] = scheduled_profile[i];
        }

        const double denom = s_max[i] - s_min[i];
        if (denom > 0.0) {
            fraction[i] = (sx[i] - s_min[i]) / denom;
        }
        if (fraction[i] < 0.0) {
            fraction[i] = 0.0;
        } else if (fraction[i] > 1.0) {
            fraction[i] = 1.0;
        }
    }

    time_t aggregated_start = fa.get_aggregated_scheduled_start_time();
    if (aggregated_start == 0) {
        aggregated_start = fa.get_aggregated_earliest();
    }

    for (auto &flexoffer : F) {
        const auto profile = flexoffer.get_profile();
        vector<double> allocations((size_t)flexoffer.get_duration(), 0.0);

        const double start_diff_sec = difftime(flexoffer.get_est(), aggregated_start);
        const int offset_hour = static_cast<int>(floor(start_diff_sec / 3600.0));

        for (int j = 0; j < flexoffer.get_duration(); j++) {
            const int idx = offset_hour + j;
            if (idx >= 0 && static_cast<size_t>(idx) < profile_size) {
                const double slice_min = profile[j].min_power;
                const double slice_max = profile[j].max_power;
                allocations[static_cast<size_t>(j)] = slice_min + (slice_max - slice_min) * fraction[idx];
            }
        }

        flexoffer.set_scheduled_allocation(allocations);
        flexoffer.set_scheduled_start_time(flexoffer.get_est());
    }
}
