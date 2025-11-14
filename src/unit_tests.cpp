#include <vector>
#include <iostream>
#include <cmath>

#include "../include/unit_test.h"
#include "../include/helperfunctions.h"
#include "../include/aggregation.h"
#include "../include/tec.h"

using namespace std;

int flexoffer_unittest();
int tec_unittest();
int nToMAggregationWithStartAlignment_unittest();
int nToMAggregationWithBalanceAlignment_unittest();
int scheduledFlexOfferDisaggregation_unittest();

int runUnitTests(){
    int parsing = 0;
    int error = 1;

    cout << "Running all tests:\n";
    if(flexoffer_unittest()){
       cout << "Flexoffer unit test passed.\n";
    } else return error;
    if(tec_unittest()){
       cout << "Tec flexoffer unit test passed.\n";
    } else return error;
    if(nToMAggregationWithStartAlignment_unittest()){
       cout << "nToMAggregation with start alignment unit test passed.\n";
    } else return error;
    if(nToMAggregationWithBalanceAlignment_unittest()){
       cout << "nToMAggregation with balance alignment unit test passed.\n";
    } else return error;
    if(scheduledFlexOfferDisaggregation_unittest()){
       cout << "Scheduled flex-offer disaggregation unit test passed.\n";
    } else return error;

    return parsing;
}

int flexoffer_unittest(){
    vector<TimeSlice> testProfile{{0, 1}};
    Flexoffer testOffer(1, 1, 2, 3, testProfile, 1);

    if(testOffer.get_offer_id() != 1){
        cout << "Flexoffer unit test has failed. Flexoffers offer id was "
             << testOffer.get_offer_id() << " but was expected to be 1.\n";
        return 0;
    }
    if(testOffer.get_est() != 1){
        cout << "Flexoffer unit test has failed. Flexoffers earliest start time was " 
             << testOffer.get_est() << " but was expected to be 1.\n";
        return 0;
    }
    if(testOffer.get_lst() != 2){
        cout << "Flexoffer unit test has failed. Flexoffers latest start time was "
             << testOffer.get_lst() << " but was expected to be 2.\n";
        return 0;
    }
    if(testOffer.get_et() != 3){
        cout << "Flexoffer unit test has failed. Flexoffers latest end time was "
             << testOffer.get_et() << " but was expected to be 3.\n";
        return 0;
    }
    if(testOffer.get_profile()[0].min_power != 0){
        cout << "Flexoffer unit test has failed. Flexoffers profile[0] was "
             << testOffer.get_profile()[0].min_power << " but was expected to be 0.\n";
        return 0;
    }
    if(testOffer.get_profile()[0].max_power != 1){
        cout << "Flexoffer unit test has failed. Flexoffers profile[1] was "
             << testOffer.get_profile()[0].max_power << " but was expected to be 1.\n";
        return 0;
    }
    if(testOffer.get_duration() != 1){
        cout << "Flexoffer unit test has failed. Flexoffers duration was "
             << testOffer.get_et() << " but was expected to be 1.\n";
        return 0;
    }

    return 1;
};

int tec_unittest(){
    vector<TimeSlice> testProfile{{0, 1}};
    Tec_flexoffer testOffer(0, 1, 1, 1, 2, 3, testProfile, 1);

    if(testOffer.get_min_overall_kw() != 0){
        cout << "TEC unit test has failed. TECs minimum kwh overall was "
             << testOffer.get_min_overall_kw() << " but was expected to be 0.\n";
        return 0;
    }
    if(testOffer.get_max_overall_kw() != 1){
        cout << "TEC unit test has failed. TECs maximum kwh overall was "
             << testOffer.get_max_overall_kw() << " but was expected to be 1.\n";
        return 0;
    }
    if(testOffer.get_offer_id() != 1){
        cout << "TEC unit test has failed. TECs offer id was " 
             << testOffer.get_offer_id() << " but was expected to be 1.\n";
        return 0;
    }
    if(testOffer.get_est() != 1){
        cout << "TEC unit test has failed. TECs earliest start time was "
             << testOffer.get_est() << " but was expected to be 1.\n";
        return 0;
    }
    if(testOffer.get_lst() != 2){
        cout << "TEC unit test has failed. TECs latest start time was "
             << testOffer.get_lst() << " but was expected to be 2.\n";
        return 0;
    }
    if(testOffer.get_et() != 3){
        cout << "TEC unit test has failed. TECs latest end time was " 
             << testOffer.get_et() << " but was expected to be 3.\n";
        return 0;
    }
    if(testOffer.get_profile()[0].min_power != 0){
        cout << "TEC unit test has failed. TECs profile[0] was " 
             << testOffer.get_profile()[0].min_power << " but was expected to be 0.\n";
        return 0;
    }
    if(testOffer.get_profile()[0].max_power != 1){
        cout << "TEC unit test has failed. TECs profile[1] was " 
             << testOffer.get_profile()[0].max_power << " but was expected to be 1.\n";
        return 0;
    }
    if(testOffer.get_duration() != 1){
        cout << "TEC unit test has failed. TECs duration was "
             << testOffer.get_et() << " but was expected to be 1.\n";
        return 0;
    }

    return 1;
};

//this is the example for Torbens Flexoffer aggregation paper from 2015 on page 2896
int nToMAggregationWithStartAlignment_unittest(){
    vector<TimeSlice> testProfile{{0, 1},{0,1}, {0,1}};
    Flexoffer testOffer1(1, 0, 4*(3600), 7*(3600), testProfile, 3);
    Flexoffer testOffer2(2, 2*(3600), 4*(3600), 6*(3600), testProfile, 3);
    vector<Flexoffer> offers{testOffer1, testOffer2};

    int est_threshold  = 4;
    int lst_threshold  = 5;
    int max_group_size = 3;
    
    vector<AggregatedFlexOffer> afos = nToMAggregation(offers, est_threshold, lst_threshold, max_group_size, Alignments::start, 0);

    if(afos.size() != 1){
        cout << "nToMAggregation_start test has failed. The size of the aggregated FO vector was "
             << afos.size() << " but was expected to be 1.\n"; 
        return 0;
    }
    if(afos[0].get_aggregated_earliest()/3600 != 0){
        cout << "nToMAggregation_start test has failed. The earliest start time of the aggregated FO vector was "
             << afos[0].get_aggregated_earliest()/3600 << " but was expected to be 0.\n"; 
        return 0;
    }
    if(afos[0].get_aggregated_latest()/3600 != 1){
        cout << "nToMAggregation_start test has failed. The latest start time of the aggregated FO vector was "
             << afos[0].get_aggregated_latest()/3600<< " but was expected to be 1.\n"; 
        return 0;
    }
    if(afos[0].get_aggregated_end_time()/3600 != 6){
        cout << "nToMAggregation_start test has failed. The latest end time of the aggregated FO vector was "
             << afos[0].get_aggregated_end_time()/3600 << " but was expected to be 6.\n"; 
        return 0;
    }
    if(afos[0].get_duration() != 5){
        cout << "nToMAggregation_start test has failed. The duration of the aggregated FO vector was "
             << afos[0].get_duration() << " but was expected to be 5.\n"; 
        return 0;
    }
    if(afos[0].get_aggregated_profile().size() != 5){
        cout << "nToMAggregation_start test has failed. The profile size of the aggregated FO vector was "
             << afos[0].get_aggregated_profile().size() << " but was expected to be 5.\n"; 
        return 0;
    }
    if(afos[0].get_aggregated_profile()[0].min_power != 0){
        cout << "nToMAggregation_start test has failed. The profile element(min power) of the aggregated FO vector was "
             << afos[0].get_aggregated_profile()[0].min_power << " but was expected to be 0.\n"; 
        return 0;
    }
    if(afos[0].get_aggregated_profile()[0].max_power != 1){
        cout << "nToMAggregation_start test has failed. The profile element(min power) of the aggregated FO vector was "
             << afos[0].get_aggregated_profile()[0].max_power << " but was expected to be 2.\n"; 
        return 0;
    }
    
    return 1;
};

int nToMAggregationWithBalanceAlignment_unittest(){
    vector<TimeSlice> testProfile{{0, 1}};
    Flexoffer testOffer1(1, 1*(3600), 2*(3600), 3*(3600), testProfile, 1);
    Flexoffer testOffer2(2, 2*(3600), 3*(3600), 4*(3600), testProfile, 1);
    vector<Flexoffer> offers{testOffer1, testOffer2};

    int est_threshold  = 2;
    int lst_threshold  = 2;
    int max_group_size = 3;
    
    vector<AggregatedFlexOffer> afos = nToMAggregation(offers, est_threshold, lst_threshold, max_group_size, Alignments::balance, 0);

    if(afos.size() != 1){
        cout << "nToMAggregation_balance test has failed. The size of the aggregated FO vector was "
             << afos.size() << " but was expected to be 1.\n";
        return 0;
    }
    if(afos[0].get_aggregated_earliest()/3600 != 2){
        cout << "nToMAggregation_balance test has failed. The earliest start time of the aggregated FO vector was "
             << afos[0].get_aggregated_earliest()/3600 << " but was expected to be 2.\n"; 
        return 0;
    }
    if(afos[0].get_aggregated_latest()/3600 != 2){
        cout << "nToMAggregation_balance test has failed. The latest start time of the aggregated FO vector was "
             << afos[0].get_aggregated_latest()/3600<< " but was expected to be 2.\n"; 
        return 0;
    }
    if(afos[0].get_aggregated_end_time()/3600 != 3){
        cout << "nToMAggregation_balance test has failed. The latest end time of the aggregated FO vector was "
             << afos[0].get_aggregated_end_time()/3600 << " but was expected to be 3.\n"; 
        return 0;
    }
    if(afos[0].get_duration() != 1){
        cout << "nToMAggregation_balance test has failed. The duration of the aggregated FO vector was "
             << afos[0].get_duration() << " but was expected to be 1.\n"; 
        return 0;
    }
    if(afos[0].get_aggregated_profile()[0].min_power != 0){
        cout << "nToMAggregation_balance test has failed. The profile element(min power) of the aggregated FO vector was "
             << afos[0].get_aggregated_profile()[0].min_power << " but was expected to be 0.\n"; 
        return 0;
    }
    if(afos[0].get_aggregated_profile()[0].max_power != 2){
        cout << "nToMAggregation_balance test has failed. The profile element(min power) of the aggregated FO vector was "
             << afos[0].get_aggregated_profile()[0].max_power << " but was expected to be 2.\n"; 
        return 0;
    }
    
    return 1;
};

int scheduledFlexOfferDisaggregation_unittest(){
    const time_t base = 1'700'000'000; // arbitrary fixed epoch base

    vector<TimeSlice> profile1{{0.0, 2.0}, {1.0, 3.0}};
    vector<TimeSlice> profile2{{0.5, 1.5}, {0.5, 2.5}};

    Flexoffer offer1(1, base, base + 1 * 3600, base + 2 * 3600, profile1, 2);
    Flexoffer offer2(2, base + 1 * 3600, base + 2 * 3600, base + 3 * 3600, profile2, 2);

    vector<Flexoffer> offers{offer1, offer2};

    AggregatedFlexOffer aggregated(100, Alignments::start, offers);
    aggregated.set_aggregated_scheduled_start_time(aggregated.get_aggregated_earliest());

    ScheduledFlexOffer scheduled(aggregated);
    scheduled.n_to_1_disaggregation(offers, aggregated);

    for (const auto &flex : offers) {
        const auto allocation = flex.get_scheduled_allocation();
        if (allocation.size() != static_cast<size_t>(flex.get_duration())) {
            cout << "Scheduled disaggregation test failed. Allocation size mismatch for offer "
                 << flex.get_offer_id() << "\n";
            return 0;
        }

        for (double value : allocation) {
            if (!std::isfinite(value)) {
                cout << "Scheduled disaggregation test failed. Non-finite allocation for offer "
                     << flex.get_offer_id() << "\n";
                return 0;
            }
        }

        if (flex.get_scheduled_start_time() != flex.get_est()) {
            cout << "Scheduled disaggregation test failed. Start time mismatch for offer "
                 << flex.get_offer_id() << "\n";
            return 0;
        }
    }

    return 1;
}
