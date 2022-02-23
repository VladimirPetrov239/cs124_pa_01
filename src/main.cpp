#include <iostream>
#include <set>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <string>
#include <math.h>
#include <bitset>
#include <map>
#include <deque>
#include <unordered_map>
#include <stack>
#include <unistd.h>
#include <queue>
#include <time.h>

using namespace std;

using ll = long long;
using ld = long double;


//working with the graph (Kruskal MST algorithm)
vector<int> parent;
vector<int> dsu_rank;

void make_set(int v) {
    dsu_rank[v] = 0;
    parent[v] = v;
}

int dsu_get(int v) {
    return (v == parent[v]) ? v : (parent[v] = dsu_get(parent[v]));
}

void dsu_unite(int u, int v) {
    u = dsu_get(u);
    v = dsu_get(v);
    if(u == v) {
        return;
    }

    if(dsu_rank[u] > dsu_rank[v]){
        swap(u, v);
    }
    parent[u] = v;
    if(dsu_rank[u] == dsu_rank[v]) {
        dsu_rank[v]++;
    }
}

ld kruskal_mst(vector<pair<ld, pair<int, int>>> g, int n) {
    parent.resize(n);
    dsu_rank.resize(n);
    for (int i = 0; i < n; ++i) {
        make_set(i);
    }

    ld cost = 0;
    int cnt = 0;
    //vector<pair<int, int>> ans; in case we wanted to store edges as well
    for(int i = 0; i < g.size(); i++) {
        pair<int, int> e = g[i].second;
        int u = e.first, v = e.second;
        if(dsu_get(u) != dsu_get(v)) {
            dsu_unite(u, v);
            //ans.push_back(make_pair(u, v));
            cnt++;
            cost += sqrt(g[i].first);
            if(cnt == n - 1) {
                // cout << "max edge = " << g[i].first <<";  ";
                return cost;
            }
        }
    }
    return cost;
}
//------------



// Generating graph for dim = 0
ld bound0(ll n){
    return 1 - std::pow(0.1, (long double) 1 / (n / 16));
}

vector<pair<ld, pair<int, int>>> generate_graph0(ll n){
    ld border = bound0(n);
    border *= 2;

    bool full_graph = n <= 1000;
    vector<pair<ld, pair<int, int>>> g;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            long double rand_weight = ((long double) rand() / (RAND_MAX));
            if(full_graph || rand_weight <= border){
                g.push_back(make_pair(rand_weight * rand_weight, make_pair(i, j)));
            }
        }
    }
    return g;
}
//------------



// Generating graph for dim >= 2
static map<int, int> max_degree{{2, 6}, {3, 12}, {4, 24}};
static map<int, ld> average_per_cell{{2, 10.0}, {3, 4.0}, {4, 1.0}};

struct cell{
    int dim, base;
    vector<int> coordinates;
    cell() : dim(0), base(0) {
        coordinates.resize(0);
    }
    cell(cell const& src) : dim(src.dim), base(src.base) {
        coordinates.resize(src.dim);
        for(int i = 0; i < src.dim; i++) {coordinates[i] = src.coordinates[i];}
    }
    cell(int dim, int base) : dim(dim), base(base) {
        coordinates.resize(dim);
    }
    int get_code(){
        int ans = 0;
        int multiply =  1;
        for (int i = 0; i < dim; ++i, multiply *= base) {
            ans += multiply * coordinates[i];
        }
        return ans;
    };
};

//return code of a cell, same as if we were dealing with base (not 10-base) counting system
int get_code(vector<int>coordinates, int dim, int base){
    int ans = 0;
    int multiply = 1;
    for (int i = 0; i < dim; ++i, multiply *= base) {
        ans += multiply * coordinates[i];
    }
    return ans;
}

//for a given point returns the cell that we will assign to this point (in short, just taking nearest 1/base fractions)
cell get_cell(vector<ld> point, int base) {
    int dim = point.size();
    cell ans = cell(dim, base);
    for(int i = 0; i < dim; i++) {
        ans.coordinates[i] = (int) (point[i] * base);
    }
    return ans;
}

// provides array of all possible shifts (each of dimension = dim): not greater than distance (in absolute value) in each direction
// e.g. get_vicinity(1, 2) = [[-1, -1], [-1, 0], [-1, 1], [0, -1], [0, 0], [0, 1], [1, -1], [1, 0], [1, 1]]
vector<vector<int>> get_vicinity(int distance, int dim) {
    vector<vector<int>> ans;
    if(dim == 1) {
        for(int j = - distance; j < distance + 1; j++){
            ans.push_back(vector<int>{j});
        }
        return ans;
    }
    vector<vector<int>> previous_level = get_vicinity(distance, dim - 1);
    vector<int> temporary;
    for(int i = 0; i < previous_level.size(); i++) {
        temporary = previous_level[i];
        for(int j = - distance; j < distance + 1; j++){
            temporary.push_back(j);
            ans.push_back(temporary);
            temporary.pop_back();
        }
    }
    return ans;
}

// same as above but leaves only "non-negative" shifts: the first non-zero coordinate change must be positive (or all shifts = 0 in case it's identical)
// e.g. get_non_negative_vicinity(1, 2) = [[0, 0], [0, 1], [1, -1], [1, 0], [1, 1]]
vector<vector<int>> get_non_negative_vicinity(int distance, int dim) {
    vector<vector<int>> old_ans = get_vicinity(distance, dim);
    vector<vector<int>> ans;
    for(int i = 0; i < old_ans.size(); i++){
        int flag_first_is_positive = 0;
        for(int k : old_ans[i]){
            if (k > 0){
                flag_first_is_positive = 1;
                break;
            }
            if(k < 0) {
                flag_first_is_positive = -1;
                break;
            }
        }
        if(flag_first_is_positive >= 0){
            ans.push_back(old_ans[i]);
        }
    }
    return ans;
}

ld square_dist (vector<ld> p1, vector<ld> p2) {
    ld sum = 0;
    for (size_t i = 0; i < p1.size(); ++i) {
        sum += (p1[i] - p2[i]) * (p1[i] - p2[i]);
    }
    return sum;
}

vector<pair<ld, pair<int, int>>> generate_graph(int n, int dim){

    vector<vector<vector<int>>> shifts;
    for (int distance = 0; distance <= 4; ++distance) {
        shifts.push_back(get_non_negative_vicinity(distance, dim));
    }

    vector<vector<ld>> points(n);
    for (int i = 0; i < n; ++i) {
        points[i].resize(dim); // points are stored here as sequence of coordinates
        for (int j = 0; j < dim; ++j) {
            points[i][j] = ((long double) rand() / (RAND_MAX));
        }
    }

    vector<pair<ld, pair<int, int>>> g;
    //generating full graph in small cases
    if(n <= 2048) {
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; ++j) {
                g.push_back(make_pair(square_dist(points[i], points[j]), make_pair(i, j)));
            }
        }
        return g;
    }

    ld avrg_per_cell = average_per_cell[dim];
    int base = (int) std::pow(n / avrg_per_cell, (long double) 1.0 / dim);
    // we will separate our box into cells, splitting evenly each side into base pieces
    // suppose we want avrg_per_cell in each cell, so we need (base ^ dim) * avrg_per_cell = n
    int cell_num = pow(base, dim);

    // we will encode all cells as vectors of length dim <-> same as using other base (than 10) for counting system
    // also, for each cell we will store numbers of points that got into this specific cell:
    vector<vector<int>> model(cell_num); // model[cell] returns number of vertexes in that cell
    // for convenience we want also to have access from the point to its cell
    vector<cell> cell_from_point(n);
    for(int i = 0; i < n; i++) {
        cell current = get_cell(points[i], base);
        model[current.get_code()].push_back(i);
        cell_from_point[i] = current;
    }

    vector<ld> distance_to_border(n); //for each point return distance to the border of the cell in which the point is situated:
    // in essence, this is the minimum (in absolute value) difference of one of the coordinates
    for (int i = 0; i < n; ++i) {
        ld cur_distance = 1;
        for (int pos = 0; pos < dim; ++pos) {
            cur_distance = min(cur_distance, abs(points[i][pos] - cell_from_point[i].coordinates[pos]));
            cur_distance = min(cur_distance, abs(points[i][pos] - cell_from_point[i].coordinates[pos] - 1));
        }
        distance_to_border[i] = cur_distance;
    }

    int limit = max_degree[dim] / 2 + 1; // now recall that in addition to looking at only close points we also want to find >= limit
    // of close neighbors to each point thus enhancing probability that we included all required edges in our graph
    // in case we cannot find many neighbor, so we deal with "outlier", we push all N - 1 edges (x, i) in the graph
    int outliers = 0; //number of "bad luck" points: where we won't find limit(dim) number of close neighbors in their vicinity

    for(int i = 0; i < n; i++) {
        // we would like to observe neighboring (with the given cell_from_point[i]) cells
        // and check whether there are enough points to find all potential edges that could appear in MST

        // "neighboring" means that we will shift the center of cell no more than +-distance in each coordinate (= 1 or 2) -- this is implemented with the shift array
        bool sufficient_data = false; // indicator whether we found >= limit neighbors

        for (int distance = 1; distance <= 4; ++distance) {

            int cnt = 0;
            vector<int> potential_edges;
            //we will store all numbers "other" of potential candidates points[other]
            //that appear in neighboring cells and close enough to our current points[i]

            //this part is the slowest
            for (int j = 0; j < shifts[distance].size(); j++) {
                vector<int> dest(dim); // coordinates of cell after shifting
                bool correct_shift = true;
                for(int pos = 0; pos < dim; pos++) {
                    dest[pos] = shifts[distance][j][pos] + cell_from_point[i].coordinates[pos];
                    if(dest[pos] < 0 || dest[pos] >= base) {
                        correct_shift = false; //cell could be shifted beyond the grid, in this case we don't add it
                        break;
                    }
                }
                if(correct_shift) {
                    int neighbor_code = get_code(dest, dim, base);
                    for (int other: model[neighbor_code]) {
                        if (other != i) {

                            ld radius = distance + distance_to_border[i] * base;
                            // precise description for this inequality can be found in the end of this function (!!!)
                            if (square_dist(points[i], points[other]) * (base * base) < radius * radius) {
                                potential_edges.push_back(other);
                                cnt++;
                            }
                        }
                    }
                }
            }
            if (cnt >= limit) { // this exactly means there are enough close points
                sufficient_data = true;
                for (int l = 0; l < potential_edges.size(); ++l) {
                    int other = potential_edges[l];
                    g.push_back(make_pair(square_dist(points[i], points[other]), make_pair(i, other)));
                }
                break; //at this point we don't need to observe greater vicinities for greater value of distance:
                //we already found limit(d) points close to our point -- main (highly probable) candidates for edges in MST
            }
        }
        if (!sufficient_data) { // we didn't find enough points so we push all possible edges
            outliers++;
            for (int other = 0; other < n; other++) {
                if (other != i) {
                    g.push_back(make_pair(square_dist(points[i], points[other]), make_pair(i, other)));
                }
            }
        }
    }
    // if we are interested, we can find number of outliers (for n = 262144 dim = 4 there are over 30 of them)
    // cout << "number of outliers = " << outliers << endl;
    return g;
}

/*
(!!!)
ld radius = distance + distance_to_border[i] * base;
 if (square_dist(points[i], points[other]) * (base * base) < radius * radius) { ... }

explanation:
 * for any point P, outside of cells shifted for no more than distance (1-4) from initial cell,
 * the distance from P to points[i] (current point) is at least (distance/base) + distance_to_border[i]
 * because there will be at least one coordinate of P that differs from the corresponding coordinate of points[i] at this value
 *[ term (distance/base) comes from the fact that P is outside of shifts,
 * plus we add distance_to_border[i] as a lower bound for how points[i] differs from its cell in that direction]
 *
 * seems redundant, but this messy variable makes a difference speeding up the RT in approx. 3 times
*/

//------------


int main(int argc, char **argv) {

    //flag 0 for initial task, flag 1 for runtime and values
    if(argc != 5) {
        cout << "incorrect number of input variables (expected 3, got: " << argc - 1 << ")" << endl;
        return 0;
    }

    if(strcmp(argv[1], "0") && strcmp(argv[1], "1") ) {
        cout << "incorrect flag" << endl;
        return 0;
    }

    int n = stoi(argv[2]), trials = stoi(argv[3]), dim = stoi(argv[4]);
    if(dim > 4 || dim == 2) {
        cout << "incorrect dimension (possible values: 0, 2, 3, 4)" << endl;
        return 0;
    }
    ld avrg = 0;
    if(!strcmp(argv[1], "0")) {
        for (int test = 0; test < trials; ++test) {
            vector<pair<ld, pair<int, int>>> g;
            if(dim > 0) {
                g = generate_graph(n, dim);
            } else {
                g = generate_graph0(n);
            }
            sort(g.begin(), g.end());
            ld ans = kruskal_mst(g, n);
            avrg += ans;
        }
        cout << avrg / trials << endl;
        return 0;
    }

    if(!strcmp(argv[1], "1")) {
        clock_t tStart = clock();
        vector<ld> values;
        map<string, double> times = {{"generate", 0}, {"sort", 0}, {"mst", 0}};
        for(int test = 0; test < trials; test++) {
            clock_t tStart1 = clock();
            vector<pair<ld, pair<int, int>>> g;
            if(dim > 0) {
                g = generate_graph(n, dim);
            } else {
                g = generate_graph0(n);
            }
            times["generate"] += (double) (clock() - tStart1) / CLOCKS_PER_SEC;
            //printf("Time for generating graph: %.2fs\n", (double) (clock() - tStart1) / CLOCKS_PER_SEC);

            clock_t tStart2 = clock();
            sort(g.begin(), g.end());
            times["sort"] += (double) (clock() - tStart2) / CLOCKS_PER_SEC;
            //printf("Time for sorting graph: %.2fs\n", (double) (clock() - tStart2) / CLOCKS_PER_SEC);

            clock_t tStart3 = clock();
            ld ans = kruskal_mst(g, n);
            times["mst"] += (double) (clock() - tStart3) / CLOCKS_PER_SEC;
            //printf("Time for find MST (Kruskal): %.2fs\n", (double) (clock() - tStart3) / CLOCKS_PER_SEC);
            values.push_back(ans);
        }

        cout << "values:" << endl;
        for (int i = 0; i < trials; ++i) {
            avrg += values[i];
            cout << values[i] << " ";
        }
        cout << endl << "Average weight = " << avrg / trials << endl;
        printf("Time statistics: \n");
        printf("Average time for generating graph: %.2fs\n", times["generate"] / trials);
        printf("Average time for sorting graph: %.2fs\n", times["sort"] / trials);
        printf("Average time for find MST (Kruskal): %.2fs\n", times["mst"] / trials);
        printf("Average time for all procedure: %.2fs\n", ((double)(clock() - tStart)/CLOCKS_PER_SEC) / trials);
        printf("Overall program time taken (for all trials): %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
        return 0;
    }
}
