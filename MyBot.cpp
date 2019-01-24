#include "hlt/game.hpp"
#include "hlt/constants.hpp"
#include "hlt/log.hpp"
#include <queue>
#include <random>
#include <ctime>
#include <map>
#include <math.h>
#include <vector>
#include <set>
#include <array>
#include <utility>
#include <cstdlib>
#include <algorithm>
#include <chrono>
#include <functional>

using namespace std;
using namespace hlt;

// & <- refers to address
// * <- points to & (address)
// int x = 5; <- variable
// int * Px = *x; <- pointer to x
// int x = 5;
// int & Rx = x; <- reference to x pointer?

#define NULL_POSITION Position(999, 999)

bool LOG_ENABLED = true;

std::string to_string(Position pos) {
	try {
		return "Position(" + to_string(pos.x) + ", " + to_string(pos.y) + ")";
	} catch (std::exception &e) {
		return e.what();
	}

}

std::string to_string(shared_ptr <Ship> ship) {
	try {
		return "\t-> Ship(id=" + std::to_string((int) ship->id) + " " + to_string(ship->position) + " cargo: " +
		       std::to_string((int) ship->halite) + ")";
	} catch (std::exception &e) {
		return e.what();
	}
}

std::string to_string(Direction move) {
	string m;
	try {
		switch (move) {
			case Direction::NORTH:
				m = "North";
				break;
			case Direction::SOUTH:
				m = "South";
				break;
			case Direction::EAST:
				m = "East";
				break;
			case Direction::WEST:
				m = "West";
				break;
			case Direction::STILL:
				m = "Still";
				break;
		}
		return m;
	} catch (std::exception &e) {
		return e.what();
	}
}

std::string to_string(MapCell *cell) {
	try {
		string m = "MapCell(" + to_string(cell->position) + " halite " + to_string(cell->halite);
		if (cell->is_occupied()) {
			m += +" ship" + to_string(cell->ship);
		}
		if (cell->has_structure()) {
			m += " structure " + to_string(cell->has_structure());
		}
		m += ")";
		return m;
	} catch (std::exception &e) {
		return e.what();
	}
}

template<typename T, typename priority_t>
struct PriorityQueue {
	typedef std::pair<priority_t, T> PQElement;
	std::priority_queue<PQElement, vector < PQElement>, greater<>> elements;

	inline bool empty() const {
		return elements.empty();
	}

	inline void put(T item, priority_t priority) {
		elements.emplace(priority, item);
	}

	T get() {
		T best_item = elements.top().second;
		elements.pop();
		return best_item;
	}
};

struct CellValue {
	MapCell *cell;
	float _distance;
	float value;
	float enemy_count;
	float halite;

	float max_enemies = (float) 9;
	float max_halite = (float) 1000;
	float max_distance;
	float min_halite;

	CellValue(MapCell *cell, int _distance, int halite) {
		this->cell = cell;
		this->_distance = _distance;
		this->halite = halite;
		this->value = compute_value();
	}

	float compute_value() {
		float value = (halite / _distance);
		//log::log("\t-> Compute Value Cell (" + to_string(this->cell->position.x) + ", " + to_string(this->cell->position.y) + ") [ Value: " + to_string(value) + " Halite: " + to_string(halite) + " Distance: " + to_string(_distance) + " ]");
		return value;
	}

	inline MapCell *get_cell() {
		return this->cell;
	}

	inline float get_dist() {
		return this->_distance;
	}

	inline float get_halite() {
		return this->halite;
	}

	inline float normalize(float v, float mx) {
		return v / mx;
	}

	inline float normalize(float v, float mx, float mn) {
		return (v - mn) / (mx - mn);
	}

};

struct BlockedEnemiesFriendlies {
	bool blocked;
	vector<MapCell *> enemies, friendlies;

	BlockedEnemiesFriendlies(bool blocked, vector<MapCell *> enemies, vector<MapCell *> friendlies) {
		this->blocked = blocked;
		this->enemies = enemies;
		this->friendlies = friendlies;
	}

	bool is_blocked(void) {
		return this->blocked;
	}

	vector<MapCell *> get_enemies(void) {
		return this->enemies;
	}

	vector<MapCell *> get_friendlies(void) {
		return this->friendlies;
	}

	inline const int get_e_cnt() {
		return (int) get_enemies().size();
	}

	inline const int get_f_cnt() {
		return (int) get_friendlies().size();
	}

	std::tuple<bool, vector < MapCell * >, vector<MapCell *>> get(void) {
		return tuple < bool, vector < MapCell * >, vector < MapCell * >> (is_blocked(), get_enemies(), get_friendlies());
	}

};

class Bot {
public:
	Bot(string attr = "", int val = 0.) {
		this->name = "New-Deathbot++";
		if (attr == "max_ships") {
			this->max_ships = val;
			this->adjust_max_ships = false;
		}
		if (attr == "min_halite") {
			this->min_halite = val;
			this->min_gain = val * .25;
			this->adjust_min_halite = false;
		}
		if (attr == "max_cost") {
			this->max_cost = val;
			this->adjust_max_cost = false;
		}
		if (attr == "max_return_halite") {
			this->max_return_halite = val;
			this->adjust_max_return_halite = false;
		}
		if (attr == "get_target_diags") {
			this->get_target_diags = val;
			this->adjust_get_target_diags = false;
		}
		if (attr == "min_drop_distance") {
			this->min_drop_distance = val;
			this->adjust_min_drop_distance = false;
		}
	}

	Game game;
	string name;
	bool adjust_max_ships = true;
    bool adjust_min_halite = true;
    bool adjust_max_cost = true;
    bool adjust_max_return_halite = true;
    bool adjust_get_target_diags = true;
    bool adjust_min_drop_distance = true;
    int max_ships;
    int min_halite;
    float min_gain;
    float max_cost;
    int max_return_halite;
    int min_drop_distance;
    bool get_target_diags;
	int last_spawn_turn = 500;
	int spawn_cnt = 0;
	int total_halite = 0;
    int low_halite = 10000;
    int max_halite = 0;
	int max_steps = 128;
	int max_drops;
	bool use_danger_zone = true;
	shared_ptr <Player> me;
	unique_ptr <GameMap> &game_map = game.game_map;
	map<int, vector < Position>> ship_paths;
	map<int, string> ship_statuses;
	map < int, MapCell * > ship_targets;
	int budget = 5000;
	bool did_spawn = false;
	bool did_drop = false;
	bool spawned_last_time = false;
	bool on_spawn = false;
    bool final_phase = false;
	bool drop_ship_set = false;
    bool save_for_drop = false;
    MapCell *max_halite_cell;
    MapCell *best_drop_target;
    chrono::time_point <chrono::steady_clock> start_time;
    vector <Direction> moves = {Direction::NORTH, Direction::SOUTH, Direction::EAST, Direction::WEST, Direction::STILL};
    vector <Position> danger_zone, base_pos;
    vector <shared_ptr<Ship>> ships, enemy_ships;
    vector<int> all_halite;
    vector<MapCell *> all_cells;
    map<int, int> halite_counts;

	void sort_ships(vector<shared_ptr<Ship>> &ships, vector <shared_ptr<Ship>> &moving_ships, vector <shared_ptr<Ship>> &returning_ships, vector<Position> &base_pos, vector<shared_ptr<Ship>> &drop_ships) {
		for (auto &ship : ships) {
			bool still = should_still(ship);
			bool movable = can_move(ship);
			ship_paths[ship->id].push_back(ship->position);
			int max_halite = 950;
			int return_halite = max_return_halite;
			if (find(base_pos.begin(), base_pos.end(), ship->position) != base_pos.end()) {
				transition(ship, "exploring");
				moving_ships.push_back(ship);
			} else if (still || !movable) {
				// Only use should_still in the early stages or if ship cant afford to move, then ships will still at the CellValue target
				transition(ship, "mining");
			} else if (ship->halite >= return_halite || ship_statuses[ship->id] == "returning") {
				returning_ships.push_back(ship);
				transition(ship, "returning");
			} else if (can_move(ship)) {
				moving_ships.push_back(ship);
				transition(ship, "exploring");
			} else {
				transition(ship, "mining");
			}
		}

	}

	void log(string s) { if (LOG_ENABLED) log::log(s); };

	void set_min_halite(int value) {
		// Keep min_gain > 15 < 30
		log("\t-> set_min_halite Value [ " + to_string(value) + " ]");
		min_halite = value;
		log("\t-> set_min_halite min_halite [ " + to_string(min_halite) + " ]");
		min_gain = (float) (min_halite * .25);
		log("\t-> set_min_halite min_gain [ " + to_string(min_gain) + " ]");
	}

	void adjust_parameters() {
		last_spawn_turn =  (int) (constants::MAX_TURNS - (int) ships.size() );
		if (adjust_max_ships) {
			max_ships = (int)game.players.size() == 2 ? 128 : 96;
		}
		if (adjust_min_halite) {
			set_min_halite(60);
		}
		if (adjust_max_cost) {
			max_cost = 200;
		}
		if (adjust_max_return_halite) {
			max_return_halite = 999;
		}
		if (adjust_get_target_diags) {
			get_target_diags = true;
		}
		if (adjust_min_drop_distance) {
			min_drop_distance = 15;
		}

	}

	void clean_statuses() {
		vector<int> ship_ids;
		for (const auto &it : me->ships) {
			auto ship = it.second;
			ship_ids.push_back(ship->id);
		}
		for (auto &it : ship_statuses) {
			int ship_id = it.first;
			string status = it.second;
			if (std::find(ship_ids.begin(), ship_ids.end(), ship_id) == ship_ids.end()) {
				ship_statuses.erase(ship_id);
			}
		}
	}

	void transition(shared_ptr <Ship> ship, string status) {
		string old_status = ship_statuses[ship->id];
		log("\t-> Transition " + to_string(ship) + " " + old_status + " - > " + status);
		ship_statuses[ship->id] = status;
	}

	BlockedEnemiesFriendlies is_blocked(Position source) {
		MapCell *source_cell = game_map->at(source);
		vector < MapCell * > cells;
		vector < MapCell * > enemies;
		vector < MapCell * > friendlies;
		int f_cnt = 0;
		int e_cnt = 0;
		int o_cnt = 0;
		bool blocked = false;
		for (auto neighbor : source.get_surrounding_cardinals()) {
			if (game_map->at(neighbor)->is_occupied()) {
				if (game_map->at(neighbor)->ship->owner == me->id) {
					f_cnt++;
					friendlies.push_back(game_map->at(neighbor));
				} else if (game_map->at(neighbor)->ship->owner != me->id) {
					e_cnt++;
					enemies.push_back(game_map->at(neighbor));
				}
			} else if (!game_map->at(neighbor)->is_occupied()) {
				o_cnt++;
			}
		}
		blocked = f_cnt + e_cnt < 4 ? false : true;
		BlockedEnemiesFriendlies bef = BlockedEnemiesFriendlies(blocked, enemies, friendlies);
		return bef;
	}

	bool is_inspired(Position source) {
		int cnt = 0;
		int fcnt = 0;
		bool ret = false;
		int sz = 4;
		vector < MapCell * > window = get_manhattan_window(source, sz);
		for (const auto &cell: window) {
			if (cnt >= 2) {
				ret = true;
			}
			if (cell->is_occupied()) {
				if (cell->ship->owner == me->id) {
					// Found at idx
					fcnt++;
					continue;
				} else {
					cnt++;
					continue;
				}
			}
		}
		return ret;
	}

	bool is_occupied_by_me(Position source) {
		return game_map->at(source)->is_occupied() && game_map->at(source)->ship->owner == me->id;
	}

	int distance(Position source, Position destination) {
		return game_map->calculate_distance(source, destination);
	}

	int winning_by() {
		if (!check_winning()) {
			return 0;
		}
		vector <shared_ptr<Player>> others = get_others();
		std::sort(others.begin(), others.end(), player_compare);
		if (!others.empty()) {
			int high_score = others[0]->halite;
			return me->halite - high_score;
		}
		return me->halite;
	}

	int losing_by() {
		if (check_winning()) {
			return -winning_by();
		}
		vector<shared_ptr<Player>> others = get_others();
		std::sort(others.begin(), others.end(), player_compare);
		if (!others.empty()) {
			int high_score = others[0]->halite;
			return high_score - me->halite;
		}
		return me->halite;
	}

	bool check_final_phase() {
		// Find the ship thats furthest from its closest base
		// Its distance to the closest base = turns to return
		// Final phase = max_turns - turns to return
		int furthest = 0;
		for (auto &ship: ships) {
			MapCell* closest = get_closest_base(ship);
			int dist = distance(ship->position, closest->position);
			if (dist > furthest) {
				furthest = dist;
			}
		}
		return game.turn_number >= constants::MAX_TURNS - (furthest + 2);
	}

	static bool player_compare(shared_ptr <Player> p1, shared_ptr <Player> p2) {
		// largest value first
		return p1->halite > p2->halite;
	}

	static bool value_sort(CellValue a, CellValue b) {
		float av = a.value;
		float bv = b.value;
		// largest value first
		return av > bv;
	}

	bool check_winning() {
		vector <shared_ptr<Player>> others = get_others();

		std::sort(others.begin(), others.end(), player_compare);
		if (!others.empty() && others[0]->halite < me->halite) {
			log("\t-> Winning TRUE");
			return true;
		} else if (!others.empty()) {
			log("\t-> Winning FALSE");
			return false;
		} else {
			log("\t-> Single Player Winning TRUE");
			return true;
		}
	}

	vector<int> get_ranks() {
		vector <shared_ptr<Player>> players = game.players;
		std::sort(players.begin(), players.end(), player_compare);
		vector<int> ranks;
		for (auto &p: players) {
			ranks.push_back(p->id);
		}
		return ranks;
	}

	void handle_drop(shared_ptr <Ship> &ship, vector <Command> *command_queue) {
		if (!game_map->at(ship)->structure) {
			int ship_halite = ship->halite;
            int cell_halite = game_map->at(ship->position)->halite;
            int total_discount = ship_halite + cell_halite;
			budget -= (constants::DROPOFF_COST - total_discount);
			command_queue->push_back(ship->make_dropoff());
			log("\t-> Dropoff Handled");
			did_drop = true;
			game_map->at(ship->position)->mark_safe();
			save_for_drop = false;
		}
	}

	void handle_spawn(vector <Command> *command_queue) {
		log("\t-> Handle Spawn");
		shared_ptr <Shipyard> shipyard = me->shipyard;
		log("\t\t-> Spawning!");
		command_queue->push_back(shipyard->spawn());
		spawn_cnt++;
		did_spawn = true;
		budget -= constants::SHIP_COST;
	}

	Position find_first_move(Position start, Position goal, map <Position, Position> came_from) {
		Position previous;
		Position current = goal;
		while (current != start && current != NULL_POSITION) {
			previous = current;
			current = came_from[current];
		}
		return previous;
	}

	float heuristic(Position current, Position start, Position goal) {
		int dx1 = current.x - goal.x;
		int dy1 = current.y - goal.y;
		int dx2 = start.x - goal.x;
		int dy2 = start.y - goal.y;
		//float cross = abs(euclidean_distance(current, goal) - euclidean_distance(start, goal));
		float cross = abs(dx1 * dy2 - dx2 * dy1);
		float heuristic = (cross * 0.01);
		return heuristic;
	}

	template<typename Position, typename cost_t>
	void astar(Position start, Position goal, map <Position, Position> &came_from, map <Position, cost_t> &cost_so_far, int &steps, int &max_steps) {
		PriorityQueue<Position, cost_t> frontier;
		vector <Position> neighbors;
		Position current;
		cost_t new_cost, priority, initial_cost;
		came_from[start] = NULL_POSITION;
		cost_so_far[start] = heuristic(start, start, goal);
		frontier.put(start, distance(start, goal));

		bool in_danger_zone = find(danger_zone.begin(), danger_zone.end(), start) != danger_zone.end();

		while (!frontier.empty()) {
			current = frontier.get();

			if (current == goal) break;

			neighbors.clear();
			get_neighbors(current, neighbors);

			for (Position next: neighbors) {
				next = game_map->normalize(next);
				//
				if (game_map->at(next)->is_occupied() || (game_map->at(next)->has_structure() && next != goal) || (find(danger_zone.begin(), danger_zone.end(), next) != danger_zone.end())) continue;

				new_cost = cost_so_far[current] + heuristic(next, start, goal);

				if (cost_so_far.find(next) == cost_so_far.end() || new_cost < cost_so_far[next]) {
					cost_so_far[next] = new_cost;
					priority = distance(next, goal);
					frontier.put(next, priority);
					came_from[next] = current;
				}
			}

			steps++;
			if (steps > max_steps) break;
		}
	}

	void get_neighbors(Position source, vector <Position> &neighbors) {
		for (auto pos: source.get_surrounding_cardinals()) neighbors.push_back(pos);
	}

	vector <Direction> get_unsafe_moves(Position source, Position destination) {
		const auto &normalized_source = game_map->normalize(source);
		const auto &normalized_destination = game_map->normalize(destination);
		const int dx = abs(normalized_source.x - normalized_destination.x);
		const int dy = abs(normalized_source.y - normalized_destination.y);
		const int wrapped_dx = game_map->width - dx;
		const int wrapped_dy = game_map->height - dy;

		vector <Direction> possible_moves;

		if (normalized_source.x < normalized_destination.x) {
			possible_moves.push_back(dx > wrapped_dx ? Direction::WEST : Direction::EAST);
		} else if (normalized_source.x > normalized_destination.x) {
			possible_moves.push_back(dx < wrapped_dx ? Direction::WEST : Direction::EAST);
		}

		if (normalized_source.y < normalized_destination.y) {
			possible_moves.push_back(dy > wrapped_dy ? Direction::NORTH : Direction::SOUTH);
		} else if (normalized_source.y > normalized_destination.y) {
			possible_moves.push_back(dy < wrapped_dy ? Direction::NORTH : Direction::SOUTH);
		}
		//log("Unsafe Moves " + to_string((int) possible_moves.size()));

		return possible_moves;
	}

	Direction reversed_move(Direction d) {
		switch (d) {
			case Direction::NORTH:
				return Direction::SOUTH;
			case Direction::SOUTH:
				return Direction::NORTH;
			case Direction::EAST:
				return Direction::WEST;
			case Direction::WEST:
				return Direction::EAST;
		}
		// Must be still
		return Direction::STILL;
	}

	Direction get_direction_to(Position source, Position target) {
		const auto &normalized_source = game_map->normalize(source);
		const auto &normalized_target = game_map->normalize(target);
		const int dx = abs(normalized_source.x - normalized_target.x);
		const int dy = abs(normalized_source.y - normalized_target.y);
		const int wrapped_dx = game_map->width - dx;
		const int wrapped_dy = game_map->height - dy;
		Direction move;
		if (normalized_source.x < normalized_target.x) {
			if (dx > wrapped_dx) {
				move = Direction::WEST;
			} else {
				move = Direction::EAST;
			}
		} else if (normalized_source.x > normalized_target.x) {
			if (dx < wrapped_dx) {
				move = Direction::WEST;
			} else {
				move = Direction::EAST;
			}
		} else if (normalized_source.y < normalized_target.y) {
			if (dy > wrapped_dy) {
				move = Direction::NORTH;
			} else {
				move = Direction::SOUTH;
			}
		} else if (normalized_source.y > normalized_target.y) {
			if (dy < wrapped_dy) {
				move = Direction::NORTH;
			} else {
				move = Direction::SOUTH;
			}
		} else {
			move = Direction::STILL;
		}

		//log("\t-> " + to_string(source) + " direction to " + to_string(target) + " = " + to_string(move));
		return move;
	}

	Direction get_safe_move(shared_ptr <Ship> &ship, Position destination) {
		Direction move;
		Position source = ship->position;
		map <Position, Position> came_from;
		map<Position, float> cost_so_far;
		Position path;
		MapCell *target;
		int dist = distance(source, destination);
		int steps = 0;
		int adjusted_max_steps = max_steps;
		int d;
		if (source == destination) {
			return Direction::STILL;
		}
		if (dist == 1 && !game_map->at(destination)->is_occupied()) {
			move = get_direction_to(ship->position, destination);
			return move;
		}

		astar(source, destination, came_from, cost_so_far, steps, max_steps);
		path = find_first_move(source, destination, came_from);

		if (path != NULL_POSITION) {
			target = game_map->at(path);

			move = get_direction_to(source, target->position);

			return move;
		}

		for (Direction dir : get_unsafe_moves(source, destination)) {
			Position n = source.directional_offset(dir);
			d = distance(n, destination);
			target = game_map->at(n);
			if (d <= dist && !target->is_occupied()) {
				return dir;
			}
		}
		return game_map->naive_navigate(ship, destination);

	}

	MapCell *get_closest_base(shared_ptr <Ship> &ship) {
		return get_closest_base(ship->position);
	}

	MapCell *get_closest_base(Position source) {
		vector < MapCell * > bases;
		MapCell *target;

		int closest = game_map->width;

		bases.push_back(game_map->at(me->shipyard));

		for (auto &it : me->dropoffs) {
			shared_ptr <Dropoff> drop = it.second;
			bases.push_back(game_map->at(drop->position));
		}
		for (auto base: bases) {
			int d = distance(base->position, source);

			if (d < closest) {
				closest = d;
				target = base;
			}
		}
		return target;
	}

	MapCell *get_closest_enemy_base(Position source) {
		vector<MapCell*> enemy_bases;
		MapCell* target;
		vector <shared_ptr<Player>> others = get_others();
		for (auto o: others) {
			enemy_bases.push_back(game_map->at(o->shipyard));
			for (auto &it: o->dropoffs) {
				shared_ptr<Dropoff> drop = it.second;
				enemy_bases.push_back(game_map->at(drop->position));
			}
		}
		int closest = game_map->width;
		for (auto base: enemy_bases) {
			int d = distance(base->position, source);
			if (d < closest) {
				closest = d;
				target = base;
			}
		}
		return target;
	}

	MapCell* get_closest_unblocked_base(Position source) {
		vector < MapCell * > bases;
        MapCell *target;

        int closest = game_map->width;

        bases.push_back(game_map->at(me->shipyard));

        for (auto &it : me->dropoffs) {
            shared_ptr <Dropoff> drop = it.second;
            bases.push_back(game_map->at(drop->position));
        }
        for (auto base: bases) {
            int d = distance(base->position, source);
            BlockedEnemiesFriendlies bef = is_blocked(base->position);
            if (bef.is_blocked()) continue;
            if (d < closest) {
                closest = d;
                target = base;
            }
        }
        return target;
	}

	MapCell* get_closest_unblocked_base(shared_ptr<Ship> &ship) {
		return get_closest_unblocked_base(ship->position);
	}

	vector <shared_ptr<Player>> get_others() {
		vector <shared_ptr<Player>> others;
		for (auto &player: game.players) {
			if (player->id != me->id) {
				others.push_back(player);
			}
		}
		return others;
	}

	int get_neighbor_count(Position source, string type = "enemy", bool include_diags = true) {
		int cnt = 0;
		vector <Position> pos;
		if (include_diags) {
			pos = {Position(source.x - 1, source.y + 1), Position(source.x + 1, source.y + 1),
		                         Position(source.x + 1, source.y - 1), Position(source.x - 1, source.y - 1)};
		}
		for (auto p : source.get_surrounding_cardinals()) {
			pos.push_back(p);
		}
		for (auto p : pos) {
			MapCell *cell = game_map->at(p);
			if (cell->is_occupied()) {
				if (cell->ship->owner == me->id && type == "friendly") {
					cnt++;
				} else if (cell->ship->owner != me->id && type == "enemy") {
					cnt++;
				} else if (type == "all") {
					cnt++;
				}
			}
		}
		return cnt;
	}

	vector<MapCell *> get_window(Position source, int size = 12) {
		vector < MapCell * > window;
		for (int y = -size; y < size + 1; y++) {
			for (int x = -size; x < size + 1; x++) {
				if (x == 0 && y == 0) continue;
				MapCell *cell = game_map->at(game_map->normalize(Position(source.x + x, source.y + y)));
				if (find(window.begin(), window.end(), cell) == window.end()) {
					window.push_back(cell);
				}
			}
		}
		return window;
	}

	vector<MapCell *> get_manhattan_window(Position source, int size = 12) {
		vector < MapCell * > window;
		for (int y = -size; y <= size; y++) {
			for (int x = -size; x <= size; x++) {
				int d = abs(x) + abs(y);
				if (d <= size && x != 0 && y != 0) {
					MapCell *cell = game_map->at(Position(source.x + x, source.y + y));
					if (find(window.begin(), window.end(), cell) == window.end()) {
						window.push_back(cell);
					}
				}
			}
		}
		return window;
	}

	vector<MapCell *> get_sorted_window(Position source, int size = 12, int max_halite = 1000) {
		vector < MapCell * > window = get_window(source, size);
		vector <CellValue> cell_values;

		for (MapCell *cell: window) {
			// current cell is always occupied by the ship...so this filters "same" cell
			if (!cell->is_occupied()) {
				// Distance away from ship... lower distance = better
				int distance_value = distance(cell->position, source);
				// Number of enemies around cell... lower enemy_count = better
				int enemy_value = get_neighbor_count(cell->position);
				// Cell halite value... higher = better
				int halite_value = cell->halite;
				CellValue d = CellValue(cell, distance_value, halite_value);
				//log::log("\t-> Cell Value [ " + to_string(d.cell) + " = " + to_string(d.value) + " ]");
				cell_values.push_back(d);
			}
		}
		if (!cell_values.empty()) {
			std::sort(cell_values.begin(), cell_values.end(), value_sort);
			vector < MapCell * > new_window;
			for (CellValue d: cell_values) {
				new_window.push_back(d.get_cell());
			}
			return new_window;
		}
		return window;
	}

	void do_move(shared_ptr <Ship> &ship, Direction move, vector <Command> *cq, bool unmark = true) {
		//log("\t-> Do move " + to_string(ship) + " ->[ " + to_string(move) + " ]<-");
		Position target = ship->position.directional_offset(move);
		if (move == Direction::STILL) {
			game_map->at(ship)->mark_unsafe(ship);
            cq->push_back(ship->stay_still());
		} else {
			if (!game_map->at(target)->is_occupied() || game_map->at(target)->ship == ship) {
				if (unmark) {
					game_map->at(ship)->mark_safe();
				}
				game_map->at(target)->mark_unsafe(ship);
				cq->push_back(ship->move(move));
			}
		}

	}

	void do_unsafe_move(shared_ptr <Ship> &ship, Direction move, vector <Command> *cq, bool unmark = true) {
		//log("\t-> Do unsafe move " + to_string(ship) + " ->[ " + to_string(move) + " ]<-");
		Position target = ship->position.directional_offset(move);
		if (move == Direction::STILL) {
			game_map->at(ship)->mark_unsafe(ship);
			cq->push_back(ship->stay_still());
		} else {
			if (unmark) {
				game_map->at(ship)->mark_safe();
			}
				game_map->at(target)->mark_unsafe(ship);
			cq->push_back(ship->move(move));
		}

	}

	bool should_still(shared_ptr <Ship> &ship, Position target) {
		int amount = ship->halite;
		float gain = still_gain(target);
		float cost = move_cost(ship->position);
		bool movable = amount >= cost;
		bool move = false;
		bool stay_still = true;
		// Structure
		if (game_map->at(target)->has_structure()) {
			return move;
		}
		// Cant
		if (!movable) {
			//log("Cant amount < cost");
            return stay_still;
        }
        // Full
		if (amount >= constants::MAX_HALITE || amount >= max_return_halite) {
			//log("Full amount > MAX");
			return move;
		// Between Full and Empty
		} else {
			int need_to_go_home = max_return_halite - amount;
			if (gain >= need_to_go_home && amount < max_return_halite) {
				return stay_still;
			}
			if (gain < min_gain) {
				return move;
			}
			return stay_still;
		}
	}

	bool should_still(shared_ptr <Ship> &ship) {
		Position target = ship->position;
		return should_still(ship, target);
	}

	bool can_move(shared_ptr <Ship> &ship) {
		Position target = ship->position;
		return can_move(ship, target);
	}

	bool can_move(shared_ptr <Ship> &ship, Position target) {
		int amount = ship->halite;
		float cost = move_cost(target);
		return amount >= cost;
	}

	float still_gain(Position source) {
		MapCell *cell = game_map->at(source);
		float gain = cell->halite * .25;
		gain += is_inspired(cell->position) ? gain * 2 : 0;
		return gain;
	}

	int turns_to_fill(shared_ptr <Ship> &ship, Position target) {
		int current = ship->halite;
		int target_halite = game_map->at(target)->halite;
		if (current >= constants::MAX_HALITE) return 0;
		int cnt = 0;
		while (current < constants::MAX_HALITE || target_halite < 1) {
			int gain = target_halite * .25;
			current += gain;
			target_halite -= gain;
			cnt++;
		}
		return cnt;
	}

	float move_cost(Position from_node) {
		return game_map->at(from_node)->halite * .1;
	}

	float move_cost(Position from_node, Position to_node) {
		float from_cost = move_cost(from_node);
		from_cost += game_map->at(to_node)->is_occupied() ? 1000 : 0;
		return from_cost;
	}

	bool can_convert(shared_ptr<Ship> &ship) {
		int total_discount = get_convert_discount(ship);
		return budget >= (constants::DROPOFF_COST - total_discount);
	}

	int get_convert_discount(shared_ptr<Ship>&ship) {
		int ship_halite = ship->halite;
		int cell_halite = game_map->at(ship->position)->halite;
		int total_discount = ship_halite + cell_halite;
		return total_discount;
	}

	float get_duration(chrono::steady_clock::time_point end_time, chrono::steady_clock::time_point start_time) {
		return (float) chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count() / (float) 1000;
	}

	void fill_ships(vector <shared_ptr<Ship>> &ships) {
		ships.clear();
		for (auto &it: me->ships) {
			auto ship = it.second;
			ships.push_back(ship);
		}
	}

	void fill_bases(vector <Position> &base_pos) {
		base_pos.clear();
		base_pos.push_back(me->shipyard->position);
		for (auto &it : me->dropoffs) {
			base_pos.push_back(it.second->position);
		}

	}

	void fill_cells(vector<MapCell *> &all_cells, vector<int> &all_halite, Position &max_halite_pos, Position &best_drop_pos, int &total_halite, int &max_halite, int &low_halite) {
		all_cells.clear();
		all_halite.clear();
		danger_zone.clear();
		halite_counts.clear();
		total_halite = 0;
		max_halite = 0;
		low_halite = 10000;
		int max_drop_halite = 0;
		int max_friendlies = 0;
		int max_friendly_halite = 0;
		for (int y = 0; y < game_map->height; y++) {
			for (int x = 0; x < game_map->width; x++) {
				Position p = Position(x, y);
				MapCell *cell = game_map->at(p);

				int halite = cell->halite;
				all_halite.push_back(halite);
				all_cells.push_back(cell);
				total_halite += halite;
				MapCell *closest_base = get_closest_base(cell->position);
                int dist = distance(cell->position, closest_base->position);

				if (halite > max_halite) {
					max_halite = halite;
					max_halite_pos = cell->position;
				}
				if (halite > max_drop_halite && dist > min_drop_distance) {
					max_drop_halite = halite;
					best_drop_pos = cell->position;
				}
				if (halite < low_halite && halite > 5) {
					low_halite = halite;
				}
				int f_neighbors = get_neighbor_count(p, "friendly", true);
				if (f_neighbors > max_friendlies) {
					max_friendlies = f_neighbors;
                    max_friendly_halite = halite;
				}
				if (cell->structure && cell->structure->owner != me->id) {
					// If cell has a shipyard or dropoff and its not mine.
					danger_zone.push_back(p);
				}
				if (cell->is_occupied() && !is_occupied_by_me(p)) {
                    danger_zone.push_back(p);
                    get_neighbors(p, danger_zone);

                }
			}
		}
		if (low_halite == 10000) {
            // No halite > 5
            low_halite = 0;
        }
	}

	MapCell *get_target(Position start, vector<MapCell *> &all_cells, vector<MapCell *> &handled_targets, int starting_halite = 0) {
		vector<CellValue> cell_values;
		MapCell* best_cell = nullptr;
		float current_gain = still_gain(start);
		int required_for_home = max_return_halite - starting_halite;
		for (auto cell: all_cells) {
			if (cell->position == start) { continue; }
			if (cell->is_occupied()) { continue; }
			bool stop = false;
			for (auto it = ship_targets.begin(); it != ship_targets.end(); ++it) {
				if (it->second == cell) stop = true;
				if (stop) break;
			}
			if (stop) continue;
			float gain = still_gain(cell->position);

            if (find(danger_zone.begin(), danger_zone.end(), cell->position) != danger_zone.end() ||
                find(handled_targets.begin(), handled_targets.end(), cell) != handled_targets.end() ||
                //find(old_pos.end()-3, old_pos.end(), cell->position) != old_pos.end() ||
                (gain < required_for_home && gain < min_gain) || gain < current_gain ) {
				continue;
			}
			float dist = distance(start, cell->position);
			float halite = gain;
			vector<Position> add_neighbors;
			get_neighbors(cell->position, add_neighbors);
			for (auto p: add_neighbors) {
				halite += still_gain(p);
			}
			CellValue d = CellValue(cell, dist, halite);
			cell_values.push_back(d);
		}
		if (!cell_values.empty()) {
			std::sort(cell_values.begin(), cell_values.end(), value_sort);
			CellValue best_cv = cell_values[0];
			best_cell = best_cv.get_cell();
			return best_cell;
		}
		return best_cell;
	}

	MapCell *get_target(shared_ptr <Ship> &ship, vector<MapCell *> &all_cells, vector<MapCell *> &handled_targets) {
        return get_target(ship->position, all_cells, handled_targets, ship->halite);
    }

    map<int, int> get_true_scores() {
        map<int, int> true_scores;
        for (auto &p: game.players) {
            int player_halite = p->halite;
            int ship_halite = 0;
            for (auto &it: p->ships) {
                auto ship = it.second;
                ship_halite += ship->halite;
            }
            true_scores[p->id] = player_halite + ship_halite;
        }
        return true_scores;
    }
	/*
    std::tuple<int, float> simulate_spawn() {
        // Need to simulate spawning on shipyard, getting a target, moving to target -> getting max_return_halite -> returning to shipyard
        // to determine if spawning is worthwhile
        int turns_left = constants::MAX_TURNS - game.turn_number;
        int current_turn = game.turn_number;
        Position current = me->shipyard->position;
        float halite_amount = 0;
        bool done = false;
        bool at_target = false;
        bool has_target = false;
        bool returning = false;
        MapCell* target;
        //log("\t-> Simulate Spawn");
        //log("\t\t-> Turns Left " + to_string(turns_left));
        vector<MapCell*> my_cells;
        for (auto cell: all_cells) {
            my_cells.push_back(new MapCell(*cell));
        }
		vector<Position> full_path;
        map<Position, float> modded_cells;
        for (auto cell: my_cells) {
            modded_cells[cell->position] = cell->halite;
        }
        while (current_turn < constants::MAX_TURNS || !done) {
            if (find(full_path.begin(), full_path.end(), current) == full_path.end()) {
                full_path.push_back(current);
            }

            if (done || current_turn >= constants::MAX_TURNS) {
                break;
            }
            if (modded_cells.find(current) == modded_cells.end()) {
                // current not in modded_cells so add it
                modded_cells[current] = game_map->at(current)->halite;
            }
            // check inspiration
            bool insp = is_inspired(current);
            // get modded_halite
            float modded_halite = modded_cells[current];
            // compute cost/gain of modded halite with inspiration
            float current_cost = modded_halite * .1;
            float current_gain = insp ? (modded_halite * .25) + (2 * (modded_halite * .25)) : (modded_halite * .25);

            if (current_gain < min_gain) {
                if (at_target) {
                    // if we are at_target -> set at_target and has_target false
                    at_target = false;
                    has_target = false;
                }
            }
            // if we are at_target and current_gain >= min_gain
            if (at_target && current_gain >= min_gain) {
                // We reached the target -> gain until filled
                if (halite_amount < max_return_halite) {
                    // if we arent full add current_gain
                    halite_amount += current_gain;
                    // deduct 25% of modded_halite from modded_halite
                    modded_cells[current] -= modded_cells[current] * .25;
                    // increment turn
                    current_turn++;
                    // start at top of loop
                    continue;
                }
            }
            if (halite_amount >= max_return_halite || current_turn + (int) full_path.size() > turns_left) {
                has_target = false;
                at_target = false;
                // we are full, set returning true
                //returning = true;
                // instead of simulating returning just add full_path.size() to current_turn to estimated return time
                current_turn += (int) full_path.size();
                for (auto p: full_path) {
                    float return_cost = modded_cells[p] * .1;
                    halite_amount -= return_cost;
                }
                done = true;
                break;
            }
			// if we cant move
	        if (current_cost > halite_amount) {
	            // add current_gain
	            halite_amount += current_gain;
	            // deduct 25% of modded_halite from modded_halite
	            modded_cells[current] -= modded_cells[current] * .25;
	            // increment turn
	            current_turn++;
	            // start at top of loop
	            continue;
	        }
	        if (!returning && !has_target) {
	            // if we arent returning get best target
	            for (auto &cell: my_cells) {
	                int modded_halite = modded_cells[cell->position];
	                cell->halite = modded_halite;
	            }
	            target = get_target(current, my_cells, halite_amount);
	            if (target == nullptr || find(full_path.begin(), full_path.end(), target->position) != full_path.end()) {
	                target = max_halite_cell;
	            }
	            has_target = true;
	            //log("\t\t-> New Target " + to_string(target->position));
	        } else if (returning && !has_target) {
	            // if we are returning get closest base to current
	            target = get_closest_base(current);
	            //log("\t\t-> Closest Base " + to_string(target->position));
	            has_target = true;
	        }
            map <Position, Position> came_from;
            map<Position, float> cost_so_far;
            int steps = 0;
            int max_steps = 32;
	        // fill came_from
	        //log("\t\t-> Target " + to_string(target->position));
	        astar(current, target->position, came_from, cost_so_far, steps, max_steps);
	        //log("\t\t-> Astar Steps " + to_string(steps));
	        // "move" current to next position
	        current = find_first_move(current, target->position, came_from);
	        // deduct cost for moving
	        halite_amount -= current_cost;
            // increment turn
            current_turn++;

	        if (current == target->position) {
	            // if we reached our target
	            at_target = true;
	            has_target = false;
	            if (returning) {
	                // if we are also returning, we are done
	                done = true;
	            }
	        }
        }
        return std::tuple<int, float>((int)current_turn - game.turn_number, halite_amount);
    }
	*/
	int play() {
		game.ready(name);
		adjust_parameters();
		int total_starting_halite = 0;
		bool adjusting = true;

		for (;;) {
			start_time = chrono::steady_clock::now();
			vector <Command> command_queue;
			chrono::time_point <chrono::steady_clock> all_cells_start, all_cells_end, ship_time_start, ship_time_end, drop_time_start, drop_time_end, returning_time_start, returning_time_end, moving_time_start, moving_time_end, danger_zone_start, danger_zone_end, dropship_time_start, dropship_time_end;
			log("\t-> Game turn " + to_string(game.turn_number));
			game.update_frame();
			me = game.me;
			budget = me->halite;
			log("\t-> Start Budget " + to_string(budget));
			// Reset bools
			final_phase = check_final_phase();
			did_spawn = false;
			did_drop = false;
			on_spawn = false;
			drop_ship_set = false;
			// Ships
			fill_ships(ships);
			clean_statuses();

			// Bases
			drop_time_start = chrono::steady_clock::now();
			fill_bases(base_pos);
			drop_time_end = chrono::steady_clock::now();
			// End Bases

			// All Cells / Danger Zone
			all_cells_start = chrono::steady_clock::now();
			danger_zone_start = chrono::steady_clock::now();
			Position max_halite_pos;
			Position best_drop_pos;
			total_halite = 0;
			low_halite = 10000;
			max_halite = 0;
			fill_cells(all_cells, all_halite, max_halite_pos, best_drop_pos, total_halite, max_halite, low_halite);
			max_halite_cell = game_map->at(max_halite_pos);
			best_drop_target = game_map->at(best_drop_pos);
			sort(all_cells.begin(), all_cells.end(), [](MapCell *a, MapCell *b) { return a->halite > b->halite; });
			sort(all_halite.begin(), all_halite.end(), [](int a, int b) { return a > b; });
			danger_zone_end = chrono::steady_clock::now();
			all_cells_end = chrono::steady_clock::now();

			int avg_halite = total_halite / (all_cells.size());

			if (total_starting_halite == 0) {
				total_starting_halite = total_halite;
			}
			float halite_perc = (float)total_halite / total_starting_halite;
            float turn_perc = (float) game.turn_number / constants::MAX_TURNS;

			int np = (int)game.players.size();
            int w = (int) game_map->width;
			set_min_halite(avg_halite*(.5 + turn_perc));


			min_drop_distance = 15 + (int)me->dropoffs.size();

			if (adjusting) {
				// turn 1 max_return halite will be first value
				// last turn max_return_halite will be first value + 100
				max_return_halite = 900 + ( turn_perc * 100. );
            }
            log("\t-> Max Return Halite [ " + to_string(max_return_halite) + " ]");

			int total_ship_count = 0;
			map<int, int> ship_counts;
			for (auto &p: game.players) {
				total_ship_count += (int)p->ships.size();
				ship_counts[p->id] = (int)p->ships.size();
			}
			int most_ships = 0;
			int most_ship_id;
			for (auto it: ship_counts) {
				auto id = it.first;
				auto count = it.second;
				if (count > most_ships) {
					most_ships = count;
					most_ship_id = id;
				}
			}


			bool have_most_ships = most_ship_id == me->id;
			int nb_ships = (int)me->ships.size();

            nb_ships = nb_ships == 0 ? 1 : nb_ships;
			bool need_to_drop = nb_ships > 14 && me->dropoffs.size() == 0 ? true : false;

            float mod = total_halite / np / nb_ships > 1000 + avg_halite ? turn_perc : turn_perc - .01;

			log("\t-> Halite % [ " + to_string(halite_perc) + " ]");
			log("\t-> Turn % [ " + to_string(turn_perc) + " ]");

			//float mod = total_halite / max_return_halite > total_ship_count ? turn_perc : turn_perc - .01;
			// If 2 player
			if (np == 2) {
				// Hard stop at 59%
                mod = mod > .58 ? .58 : mod;
                if (w == 64) {
                    mod = mod > .54 ? .54 : mod;
                }
                if (w == 56) {
                    mod = mod > .53 ? .53 : mod;
                }
                if (w == 48) {
                    mod = mod > .52 ? .52 : mod;
                }

                if (w == 40) {
                    mod = mod > .47 ? .47 : mod;
                }
                if (w == 32) {
                    mod = mod > .45 ? .45 : mod;
                }

			}
			// If 4 Player
			if (np == 4) {
				// Hard Stop at 48%
				mod = mod > .48 ? .48 : mod;
				if (w == 40 ) {
					mod = mod > .47 ? .47 : mod;
				}
				if (w == 32) {
					mod = mod > .46 ? .46 : mod;
				}
			}

			if (!have_most_ships && np == 2) {
				mod += .01;
			}
			//bool check_perc = true;
			/*
			if (game.turn_number > 10 && !final_phase && (!have_most_ships || !check_winning()))  {
				auto sim = simulate_spawn();
				// sim = (last_turn, halite_amount)
				int turns_left = constants::MAX_TURNS - game.turn_number;
				int turns_to_deliver = std::get<0>(sim);
				float delivery_amount = std::get<1>(sim);
				if (turns_left < turns_to_deliver) {
					// if it will take more turns to deliver than we have left dont spawn
					last_spawn_turn = game.turn_number;
					check_perc = false;
				} else if (delivery_amount > losing_by() || delivery_amount > avg_halite){
					last_spawn_turn = game.turn_number + turns_to_deliver;
					check_perc = false;
				}
			}
			*/
			int turns_left = constants::MAX_TURNS - game.turn_number;
			// Hard stop at < 35% halite remaining || mod
			if (halite_perc <= 0.35 || turns_left <= 100 || turn_perc > mod) {
				last_spawn_turn = game.turn_number;
			}

			log("\t-> Last Spawn Turn [ " + to_string(last_spawn_turn) + " ]");

			switch ((int)game.players.size()) {
				case 2:
					switch (game_map->width) {
                        case 32:
                            max_drops = 2;
                            break;
                        case 40:
                            max_drops = 3;
                            break;
                        case 48:
                            max_drops = 3;
                            break;
                        case 56:
                            max_drops = 4;
                            break;
                        case 64:
                            max_drops = 8;
                            break;
                    }
                    break;
                case 4:
                    switch (game_map->width) {
                        case 32:
                            max_drops = 2;
                            break;
                        case 40:
                            max_drops = 3;
                            break;
                        case 48:
                            max_drops = 3;
                            break;
                        case 56:
                            max_drops = 4;
                            break;
                        case 64:
                            max_drops = 8;
                            break;
                    }
                    break;
			}
			if (!final_phase) {
				vector <shared_ptr<Ship>> moving_ships, returning_ships, moved_already, drop_ships;
				vector < MapCell * > handled_targets;

				ship_time_start = chrono::steady_clock::now();
				sort_ships(ships, moving_ships, returning_ships, base_pos, drop_ships);
				ship_time_end = chrono::steady_clock::now();

				returning_time_start = chrono::steady_clock::now();
				// fix ships blocked on base

				for (auto &ship: moving_ships) {
				    if (find(moved_already.begin(), moved_already.end(), ship) != moved_already.end()) {
				        // skip ships that already moved
				        continue;
				    }
				    if (find(base_pos.begin(), base_pos.end(), ship->position) != base_pos.end()) {
				        // Ship on dropoff or shipyard
				        BlockedEnemiesFriendlies ship_bef = is_blocked(ship->position);
                        // moving ship currently blocked and on base_pos
                        bool done = false;
                        if (ship_bef.is_blocked()){
                            for (auto f_cell: ship_bef.get_friendlies()) {
                                auto other_ship = f_cell->ship;
                                if (done) continue;
                                if (find(moved_already.begin(), moved_already.end(), other_ship) == moved_already.end() && can_move(other_ship)
                                    && find(returning_ships.begin(), returning_ships.end(), other_ship) != returning_ships.end()) {
                                    // other_ship hasn't moved yet and is returning so swap
                                    Direction move = get_direction_to(ship->position, f_cell->position);
                                    Direction other_move = reversed_move(move);
                                    do_unsafe_move(ship, move, &command_queue, false);
                                    do_unsafe_move(other_ship, other_move, &command_queue, false);
                                    moved_already.push_back(ship);
                                    moved_already.push_back(other_ship);
                                    done = true;
                                }

                            }
                        }
				    }
				}
				/*
				if (!drop_ships.empty()) {
					log("\t-> Drop Ships size " + to_string((int)drop_ships.size()));
					shared_ptr<Ship> dropper = drop_ships[0];

					if (can_convert(dropper)) {
                        handle_drop(dropper, &command_queue);
					} else {
						Direction move = get_direction_to(dropper->position, best_drop_target->position);
                        do_move(dropper, move, &command_queue);
                        moved_already.push_back(dropper);
					}
				}
				*/
				if (!returning_ships.empty()) {
					log("\t-> Returning ships size " + to_string((int) returning_ships.size()));

					for (auto &ship: returning_ships) {
						if (find(moved_already.begin(), moved_already.end(), ship) != moved_already.end()) {
                            // skip ships that already moved
                            continue;
                        }
						Direction move;
						Position target_pos;
						Position ship_pos = ship->position;
						// closest base blocked or unblocked
						MapCell *closest = get_closest_base(ship);
						// closest dist
						int closest_dist = distance(ship_pos, closest->position);

						// if closest is only 1 cell away
						if (closest_dist == 1) {
							if (closest->is_occupied()){
							    if (!is_occupied_by_me(closest->position)) {
									move = get_direction_to(ship->position, closest->position);
									do_unsafe_move(ship, move, &command_queue);
									moved_already.push_back(ship);
									if (closest->position == me->shipyard->position) {
										on_spawn = true;
									}

									continue;
								}
								if (is_occupied_by_me(closest->position)) {
									// Occupied by me -> try a swap
									auto docked_ship = closest->ship;
									if (find(moved_already.begin(), moved_already.end(), docked_ship) == moved_already.end()) {
										// ship hasnt moved
										MapCell* ship_cell = game_map->at(ship);
										move = get_direction_to(ship->position, docked_ship->position);
										Direction docked_move = reversed_move(move);
										do_unsafe_move(ship, move, &command_queue, false);
										do_unsafe_move(docked_ship, docked_move, &command_queue, false);
										moved_already.push_back(ship);
										moved_already.push_back(docked_ship);
										if (closest->position == me->shipyard->position) {
                                            on_spawn = true;
                                        }
										continue;
									}
								}
							} else {
							// Not occupied
								move = get_direction_to(ship->position, closest->position);
								do_move(ship, move, &command_queue);
								moved_already.push_back(ship);
								if (closest->position == me->shipyard->position) {
									on_spawn = true;
								}
								continue;
							}

						}
						// if closest is only 2 cell away
						if (closest_dist == 2) {
							// To prevent someone blocking all neighbors
							BlockedEnemiesFriendlies bef = is_blocked(closest->position);
							if ((bef.is_blocked() && bef.get_e_cnt() > 0)) {
								// Base is blocked by at least 1 enemy
								vector<MapCell*> enemies = bef.get_enemies();
								move = get_direction_to(ship->position, enemies[0]->position);
								do_unsafe_move(ship, move, &command_queue);
								moved_already.push_back(ship);
								continue;
							}
						}
						if (closest_dist > min_drop_distance && (int) me->dropoffs.size() < max_drops && (int) me->ships.size() > 12) {
							if (need_to_drop && can_convert(ship)) {
								handle_drop(ship, &command_queue);
								continue;
							}
                            int enemy_count = get_neighbor_count(ship->position, "enemy", true);
                            if (me->dropoffs.size() == 0 && can_convert(ship)) {
                                handle_drop(ship, &command_queue);
                                continue;
                            }
                            if ((enemy_count > 2 || get_convert_discount(ship) > avg_halite + 1000) && can_convert(ship)) {
                                handle_drop(ship, &command_queue);
                                continue;
                            }
                        }
						move = get_safe_move(ship, closest->position);
						do_move(ship, move, &command_queue);
						moved_already.push_back(ship);
					}
				}

				returning_time_end = chrono::steady_clock::now();
				moving_time_start = chrono::steady_clock::now();

				if (!moving_ships.empty()) {
					log("\t-> Moving ships size " + to_string(moving_ships.size()));

					for (auto &ship : moving_ships) {
						auto check_time = chrono::steady_clock::now();
                        float diff = (float) 2.0 - (float) get_duration(check_time, start_time);
                        log::log("\t-> Time Left " + to_string(diff));
                        if (diff < 0.1) {
                            log(to_string(ship) + " Out of time!");
                            continue;
                        }
						if (find(moved_already.begin(), moved_already.end(), ship) != moved_already.end()) {
                            // skip ships that already moved
                            continue;
                        }

						Position ship_pos = ship->position;
						Direction move;
						MapCell * target;

						if (ship_targets.find(ship->id) != ship_targets.end()) {
							// Existing target
							target = ship_targets[ship->id];
							float target_gain = still_gain(target->position);
							BlockedEnemiesFriendlies bef = is_blocked(target->position);
							if (target_gain < min_gain ||
								target_gain < still_gain(ship->position) ||
								(find(danger_zone.begin(), danger_zone.end(), target->position) != danger_zone.end()) ||
								bef.is_blocked()){
								ship_targets.erase(ship->id);
							} else {
								move = get_safe_move(ship, target->position);
								do_move(ship, move, &command_queue);
								handled_targets.push_back(target);
								moved_already.push_back(ship);
								continue;
							}
						}

						target = get_target(ship, all_cells, handled_targets);

                        if (target != nullptr) {
                            handled_targets.push_back(target);
                            ship_targets[ship->id] = target;
                            move = get_safe_move(ship, target->position);
                            do_move(ship, move, &command_queue);
                            moved_already.push_back(ship);
                            continue;
                        }

						log("\t-> No Targets");
						target = max_halite_cell;
						move = get_safe_move(ship, target->position);
						do_move(ship, move, &command_queue);
						moved_already.push_back(ship);
						continue;
					}
				}

				moving_time_end = chrono::steady_clock::now();

			} else {
				log("\t-> Final Phase [ True ]");
				for (auto &it : me->ships) {
					Direction move;
					shared_ptr <Ship> ship = it.second;
					MapCell *target = get_closest_base(ship);
					int dist = distance(ship->position, target->position);
					if (dist == 0) {
						if (dist == 0) {
                            // already on base - check for enemies
                            log("\t-> Ship " + to_string(ship->id) + " dist == 0 final phase");
                            BlockedEnemiesFriendlies bef = is_blocked(ship->position);
                            int max_enemy_amt = 0;
                            MapCell* max_enemy_cell;
                            for (auto p: bef.get_enemies()) {
                                if (p->ship->halite > max_enemy_amt) {
                                    max_enemy_amt = p->ship->halite;
                                    max_enemy_cell = p;
                                }
                            }
                            if (max_enemy_amt > 0) {
                                move = get_direction_to(ship->position, max_enemy_cell->position);
                                do_unsafe_move(ship, move, &command_queue);
                                continue;
                            }
                        }
					}
                    if (dist >= 2) {
                        move = get_safe_move(ship, target->position);
                        do_move(ship, move, &command_queue);
                        continue;
                    }
                    if (dist == 1) {
                        move = get_direction_to(ship->position, target->position);
                        do_unsafe_move(ship, move, &command_queue);
                        continue;
                    }
				}
			}
			// End Moving Ships
			int wby = winning_by();
			if (game.turn_number <= 1) {
				log("\t-> Game.turn_number <= 1");
				log("\t-> Max Ships [ " + to_string(max_ships) + " ]");
				handle_spawn(&command_queue);
			} else if (budget >= constants::SHIP_COST && !did_spawn) {
				if (!final_phase) {
					bool can_spawn = true;
					BlockedEnemiesFriendlies bef = is_blocked(me->shipyard->position);
					bool blocked = bef.is_blocked();
					bool occ = game_map->at(me->shipyard)->is_occupied();
					bool is_occ_by_me = is_occupied_by_me(me->shipyard->position);
					if ((occ && is_occ_by_me)) {
						can_spawn = false;
					}
					if (can_spawn) {
						if (spawn_cnt < max_ships && game.turn_number < last_spawn_turn) {
                            handle_spawn(&command_queue);
                        }
					}
				}
			}

			spawned_last_time = did_spawn;
			string spawned = spawned_last_time ? "True" : "False";
			string finals =  final_phase ? "True" : "False";

			auto end_time = chrono::steady_clock::now();
			log("\t-> Turn time [ " + to_string(get_duration(end_time, start_time)) + " ]");
			log("\t-> Ship time [ " + to_string(get_duration(ship_time_end, ship_time_start)) + " ]");
			log("\t-> Moving Ship time [ " + to_string(get_duration(moving_time_end, moving_time_start)) + " ]");
			log("\t-> Drop Ship time [ " + to_string(get_duration(drop_time_end, drop_time_start)) + " ]");
			log("\t-> Returning Ship time [ " + to_string(get_duration(returning_time_end, returning_time_start)) + " ]");
			log("\t-> Drop Time [ " + to_string(get_duration(drop_time_end, drop_time_start)) + " ]");
			log("\t-> Danger Zone Time [ " + to_string(get_duration(danger_zone_end, danger_zone_start)) + " ]");
			log("\t-> All Cell Time [ " + to_string(get_duration(all_cells_end, all_cells_start)) + " ]");
			log("\t-> # Ships [ " + to_string(ships.size()) + " ]");
			log("\t-> Spawned [ " + spawned + " ]");
			log("\t-> Winning By [ " + to_string(wby) + " ]");
			log("\t-> Final Phase [ " + finals + " ]");
			log("\t-> End Budget [ " + to_string(budget) + " ]");
			log("\t-> Min Gain [ " + to_string(min_gain) + " ]");
			log("\t-> Avg Halite [ " + to_string(avg_halite) + " ]");
			log("\t-> Last Spawn Turn [ " + to_string(last_spawn_turn) + " ]");

			if (!game.end_turn(command_queue)) {
				break;
			}

		};
		return 0;
	}

};

int main(int argc, char *argv[]) {
	string attr_str;
	int val;
	if (argc > 1) {
		attr_str = argv[1];
		val = atoi(argv[2]);
		Bot bot{attr_str, val};
		try {
			bot.play();
		} catch (exception &e) {
			log::log("\t-> Error ");
			log::log(e.what());
		}
	} else {
		Bot bot{};
		try {
			bot.play();
		} catch (exception &e) {
			log::log("\t-> Error ");
			log::log(e.what());
		}
	}
	return 0;
}
