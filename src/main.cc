//#define _POSIX_SOURCE
#include <boost/regex.hpp>
#include <boost/program_options.hpp>
#include <ctime>
#include <experimental/filesystem>
#include <fstream>
#include <exception>
#include <limits>
#include <signal.h>
#include <unistd.h>
#include <stdio.h> 
#include <sys/ipc.h> 
#include <sys/shm.h> 
#include <sys/wait.h>
#include <chrono>

#ifdef RELATION
#include "../include/correlated_sampling.h"
#include "../include/bound_sketch.h"
#else
#ifdef INTERSECTION
#include "../include/pattern_mining.h"
#include "../include/alley.h"
#else
#include "../include/cset.h"
#include "../include/impr.h"
#include "../include/sumrdf.h"
#include "../include/wander_join.h"
#include "../include/index_sample.h"
#include "../include/jsub.h"
#endif
#endif
#include "../include/memory.h"

namespace po = boost::program_options;
typedef std::numeric_limits< double > dbl;
typedef std::chrono::high_resolution_clock Clock;

struct QueryResult {
  double est;
  double time;
  int m_est;
};

struct ThreadInput {
    DataGraph* g;
    QueryGraph* q;
    double p;
    Estimator* estimator;
    int seed;
    QueryResult result;

    pthread_cond_t cond;
    pthread_mutex_t mutex;
};

void* t_function(void* data) {
    ThreadInput* input = (ThreadInput*) data;
#ifndef GGDB
    int oldtype;
#endif

    srand(input->seed);
    auto chkpt = Clock::now();
    input->result.est = input->estimator->Run(input->g, input->q, input->p);
    auto elapsed_nanoseconds = Clock::now() - chkpt;
#ifdef PRINT_STAT
	input->estimator->PrintStat();
#endif
    input->result.time = (double) elapsed_nanoseconds.count() / 1000000;
    input->result.m_est = std::max(input->result.m_est, getValueOfPhysicalMemoryUsage());

#ifndef GGDB
    const int lock_rv = pthread_mutex_lock(&(input->mutex));
    const int signal_rv = pthread_cond_signal(&(input->cond));
    const int unlock_rv = pthread_mutex_unlock(&(input->mutex));
#endif

    return NULL;
}

int main(int argc, char** argv) {

    po::options_description desc("gCare Framework");
    desc.add_options()
        ("help,h", "Display help message")
        ("query,q", "query mode")
        ("build,b", "build mode")
        ("method,m", po::value<std::string>(), "estimator method")
        ("input,i", po::value<std::string>(), "input file (in build mode: text data graph, in query mode: text query graph)")
        ("output,o", po::value<std::string>(), "output directory in query mode")
        ("data,d", po::value<std::string>(), "binary datafile")
        ("ratio,p", po::value<string>()->default_value("0.03"), "sampling ratio")
        ("iteration,n", po::value<int>()->default_value(30), "iterations per query")
        ("seed,s", po::value<int>()->default_value(0), "random seed");
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);

    if (vm.count("help") || !vm.count("input") || !vm.count("data")) {
        cout << desc;
        return -1;
    }

    if (!vm.count("query") && !vm.count("build")) {
        cout << "mode is not specified" << endl;
        cout << desc;
        return -1;
    }

    if (vm.count("query") && vm.count("build")) {
        cout << "only one mode can be set" << endl;
        cout << desc;
        return -1;
    }

    const string input_str = vm["input"].as<string>();
    const string data_str  = vm["data"].as<string>();
    double p      = stod(vm["ratio"].as<string>()); 
    int seed      = vm["seed"].as<int>();

    string method = vm["method"].as<string>(); 
    Estimator* estimator = nullptr;
    string summary_str = data_str + string(".") + method;

#ifdef RELATION
    if (method == string("cs")) {
        summary_str = summary_str + ".p" + vm["ratio"].as<string>();
        estimator = new CorrelatedSampling;
    } else if (method == string("bsk")) {
        summary_str = summary_str + ".b" + string(getenv("GCARE_BSK_BUDGET"));
        estimator = new BoundSketch;
        p = std::stod(string(getenv("GCARE_BSK_BUDGET")));
    }
#else
#ifdef INTERSECTION
    if (method == string("alley")) {
        summary_str = summary_str + ".p" + vm["ratio"].as<string>();
        estimator = new Alley(false);
    } else if (method == string("alleyTPI")) {
        summary_str = summary_str + ".L" + string(getenv("GCARE_ALLEY_TPI_MAXL")) + ".g" + string(getenv("GCARE_ALLEY_TPI_NUM_GROUP")) + ".t" + string(getenv("GCARE_ALLEY_TPI_FAIL_THRESHOLD"));
        estimator = new Alley(true);
    }
#else
    if (method == string("cset")) {
        estimator = new CharacteristicSets;
    } else if (method == string("impr")) {
        summary_str = summary_str + ".p" + vm["ratio"].as<string>();
        estimator = new Impr;
    } else if (method == string("sumrdf")) {
        summary_str = summary_str + ".p" + vm["ratio"].as<string>() + ".t" + string(getenv("GCARE_SUMRDF_THRESHOLD"));
        estimator = new SumRDF;
    } else if (method == string("wj")) {
        summary_str = summary_str + ".p" + vm["ratio"].as<string>();
        estimator = new WanderJoin;
    } else if (method == string("jsub")) {
        summary_str = summary_str + ".p" + vm["ratio"].as<string>();
        estimator = new JSUB;
    } else if (method == string("ibjs")) {
        summary_str = summary_str + ".p" + vm["ratio"].as<string>();
        estimator = new IndexSample;
    }
#endif
#endif

    summary_str = summary_str + ".s" + to_string(seed);

  std::cout << "summary: " << summary_str << "\n";
  string output_str  = vm["output"].as<string>();

  if (vm.count("build")) {
    //build mode
    DataGraph g;
    if (!g.HasBinary(data_str.c_str())) {
      std::cout << "There is no binary\n";
      g.ReadText(input_str.c_str());
      g.MakeBinary();
      g.WriteBinary(data_str.c_str());
      g.ClearRawData();
    }
    g.ReadBinary(data_str.c_str());

    srand(seed);
    auto chkpt = Clock::now();
    estimator->Summarize(g, summary_str.c_str(), p); 
    auto elapsed_nanoseconds = Clock::now() - chkpt;
    double summary_build_time = elapsed_nanoseconds.count() / 1000000;
    std::fstream fout;
    string output_fn = output_str;
    fout.open(output_fn.c_str(), std::fstream::out);
    fout << summary_build_time << endl;
    fout.close();
  } else {
      int tics_per_second = sysconf(_SC_CLK_TCK);
    //query mode
    std::cout << "Query Mode\n";
    {
        DataGraph g;
        int chkpt = getValueOfPhysicalMemoryUsage();
        g.ReadBinary(data_str.c_str());
        estimator->ReadSummary(summary_str.c_str());
        int num_iter  = vm["iteration"].as<int>();
        namespace fs = std::experimental::filesystem;
        std::fstream fout;
        string output_fn = output_str; 
        fout.open(output_fn.c_str(), std::fstream::out);
        string err_fn = output_str + string(".err");
        std::fstream err_fout;
        err_fout.open(err_fn.c_str(), std::fstream::out);

        //Add for timeout
        vector<string> timeout_group; // EXIT FOR TIMEOUT GROUP
        timeout_group.clear();
        string timeout_fn = output_str + string(".timeout");
        std::fstream timeout_fout;
        timeout_fout.open(timeout_fn.c_str(), std::fstream::out);
        string rootfiledir = input_str;
        for (auto& dir_entry : fs::recursive_directory_iterator(input_str.c_str())) {
            std::cout << "Start " << dir_entry.path().string() << "\n";
            if (dir_entry.path().string().find_last_of(".txt") + 1 != dir_entry.path().string().length()) continue;

            string group_name("NOGROUP");

            string cur_file_path=dir_entry.path().string();
            string cur_file_name=cur_file_path.substr( cur_file_path.find_last_of("/")+1 );
            string cur_file_prefix1=cur_file_path.substr(0, cur_file_path.find_last_of("/")+1 );
            string cur_file_prefix2=cur_file_path.substr(0, cur_file_path.find_last_of("/") );
            if(cur_file_prefix1 != rootfiledir && cur_file_prefix2 != rootfiledir) {
                // Assumption: If a path has a subdir, queries are grouped and the nae of the subdir is the name of the group
                group_name=cur_file_prefix2;
                if(cur_file_prefix2.find_last_of("/") != string::npos)
                    group_name=cur_file_prefix2.substr( cur_file_prefix2.find_last_of("/") + 1 );
            }
            if(std::find(timeout_group.begin(), timeout_group.end(), group_name) != timeout_group.end()) {
                std::cout << "This group " << group_name << "is already in timeout group; then Skip for " << dir_entry.path().string() << "\n";
                continue;
            }

            std::cout << "Estimator for " << dir_entry.path().string() << "\n";
#ifndef GGDB
            try {
#endif
                QueryGraph q;
                q.ReadText(dir_entry.path().string().c_str());
                vector<double> est_vec;
                double avg_est = 0.0, avg_time = 0.0;
                int num_est = 0;
                for (int i = 0; i < num_iter; i++) {
                    ThreadInput thread_input;
                    thread_input.g = &g;
                    thread_input.q = &q;
                    thread_input.p = p;
                    thread_input.estimator = estimator;
                    thread_input.seed = seed + i;
#ifndef GGDB
                    pthread_cond_init(&thread_input.cond, NULL);
                    pthread_mutex_init(&thread_input.mutex, NULL);


                    pthread_t p_thread;
                    int status;
                    int thr_id = pthread_create(&p_thread, NULL, t_function, (void*)&thread_input); 

                    struct timespec max_wait = {0, 0};
                    clock_gettime(CLOCK_REALTIME, &max_wait);
                    max_wait.tv_sec += 60; //wait for 60 seconds 

                    pthread_mutex_lock(&thread_input.mutex);
                    const int timed_wait_rv = pthread_cond_timedwait(&thread_input.cond, &thread_input.mutex, &max_wait);
										pthread_mutex_unlock(&thread_input.mutex);
                    if (timed_wait_rv == ETIMEDOUT) {
                        std::cout << "timeout\n";
                        pthread_cancel(p_thread);
                        usleep(1000000); // sleep 1 second
                        pthread_join(p_thread, (void**) &status);

                        //add current group of queries into a list excluded
                        std::cout << "throw exception\n";
                        throw Estimator::ErrCode::TIMEOUT;
                    }
                    else if (timed_wait_rv) {
                        assert(false);
                    }
                    else {
                        pthread_join(p_thread, (void**) &status);
                        if (thread_input.result.est > -1e9) {
                            est_vec.push_back(thread_input.result.est);
                            avg_time += thread_input.result.time;
                        }
                    }
#else
                    t_function((void*)&thread_input);
                    if (thread_input.result.est > -1e9) {
                        est_vec.push_back(thread_input.result.est);
                        avg_time += thread_input.result.time;
                    }
#endif
                }
                for (double est : est_vec) avg_est += est;
                avg_est /= est_vec.size();
                avg_time /= est_vec.size();
                double var = 0.0;
                for (double est : est_vec) var += (est - avg_est) * (est - avg_est);
                var /= est_vec.size();
                int precision = std::numeric_limits<double>::max_digits10;
                cout.precision(dbl::max_digits10);
                fout.precision(dbl::max_digits10);
                fout << dir_entry.path().string() << " " << avg_est << " " << avg_time << " " << var << endl;
                cout << dir_entry.path().string() << " " << avg_est << " " << avg_time << " " << var << endl;
                fout << est_vec.size();
                cout << est_vec.size();
                for (double est : est_vec) fout << " " << est;
                for (double est : est_vec) cout << " " << est;
                fout << endl;
                cout << endl;
#ifndef GGDB
            } catch (Estimator::ErrCode e) {
                err_fout << dir_entry.path().string() << " error with code " << e << "\n";
                if(e == Estimator::ErrCode::TIMEOUT) {
                    if (group_name != string("NOGROUP")) {
                        timeout_fout << group_name << endl;
                        timeout_group.push_back(group_name);
                    }
                }
            } catch (int e) {
                err_fout << dir_entry.path().string() << " error with the signal " << e << "\n";
            }
#endif
        }
        fout.close();
        err_fout.close();
        timeout_fout.close();
        estimator->End();
    }
  }
  return 0;
}
