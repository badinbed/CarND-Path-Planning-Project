#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;
using namespace Eigen;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
        return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}


int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

        double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

double mphToMs(double mph) {
    return mph/2.24;
}

double msToMph(double ms) {
    return ms*2.24;
}

void caluclateEffectiveLaneSpeeds(vector<double>& out_lane_speeds, const vector<vector<double>>& sensor_fusion, double ref_s, double seconds_look_ahead, double safety_distance) {
    for(int object_id = 0; object_id < sensor_fusion.size(); ++object_id) {
        auto& object = sensor_fusion[object_id];
        double object_d = object[6];
        int object_lane = static_cast<int>((object_d - 2.0) / 4.0 + 0.5);

        double object_s = object[5];
        double object_vx = object[3];
        double object_vy = object[4];
        double object_speed = sqrt(object_vx*object_vx + object_vy*object_vy);
        object_s += seconds_look_ahead * object_speed;

        if(object_s > ref_s - safety_distance) {
            out_lane_speeds[object_lane] = min(msToMph(object_speed), out_lane_speeds[object_lane]);
        }
    }
}

void calculateEfficiencyCost(vector<double>& out_cost, const vector<double>& lane_speeds, int ref_lane, double target_speed) {
   int num_lanes = lane_speeds.size();
    out_cost.resize(num_lanes);
    double max_lane_speed = 0;
    int fastest_lane = numeric_limits<int>::max();
    for(int i = 0; i < num_lanes; ++i) {
        if(lane_speeds[i] > max_lane_speed || (max_lane_speed == lane_speeds[i] && abs(fastest_lane - ref_lane) > abs(i - ref_lane))) {
            max_lane_speed = lane_speeds[i];
            fastest_lane = i;
        }
        double lane_change_cost = i == ref_lane ? 0.0 : 2.0;
        out_cost[i] = (lane_change_cost + target_speed - lane_speeds[i]) / (target_speed+2.0);
    }

    // a lane might be individually insufficient but bring us closer to a more efficient lane
    // we increase the cost depending on the distance of the fastest lane to promote intermediate lane changes
    if(fastest_lane != ref_lane) {
        double potential_gain = max_lane_speed - lane_speeds[ref_lane];
        double d_weight = pow(num_lanes-1, 2);
        for(int i = 0; i < num_lanes; ++i) {
            double d = abs(fastest_lane - i) / d_weight;
            double d_cost = (d * potential_gain) / (d_weight * target_speed);
            out_cost[i] = (out_cost[i] + d_cost);
        }
    }
}

void calculateSafetyCost(vector<double>& out_safety_cost, const vector<vector<double>>& sensor_fusion, int ref_lane, double ref_s, double seconds_look_ahead, double safety_distance) {
    out_safety_cost.resize(3, 0.0);
    for(int i = 0; i < 3; ++i) {
        if(abs(i - ref_lane) > 1)
           out_safety_cost[i] = 1.0;
    }
    safety_distance *= 0.7;
    for(int object_id = 0; object_id < sensor_fusion.size(); ++object_id) {
        auto& object = sensor_fusion[object_id];
        double object_d = object[6];
        int object_lane = static_cast<int>((object_d - 2.0) / 4.0 + 0.5);
        if(abs(object_lane - ref_lane) == 1) {
            double object_s = object[5];
            double object_vx = object[3];
            double object_vy = object[4];
            double object_speed = sqrt(object_vx*object_vx + object_vy*object_vy);
            object_s += seconds_look_ahead * object_speed;

            if(object_s + safety_distance > ref_s && object_s - safety_distance < ref_s) {
                out_safety_cost[object_lane] = 1.0;
            }
        }
    }
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_x;
  vector<double> map_y;
  vector<double> map_s;
  vector<double> map_dx;
  vector<double> map_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
        istringstream iss(line);
        double x;
        double y;
        float s;
        float d_x;
        float d_y;
        iss >> x;
        iss >> y;
        iss >> s;
        iss >> d_x;
        iss >> d_y;
        map_x.push_back(x);
        map_y.push_back(y);
        map_s.push_back(s);
        map_dx.push_back(d_x);
        map_dy.push_back(d_y);
  }

  int ref_lane = 1;
  double ref_vel = 0.0; //mph
  double target_vel = 49.5; //mph
  int horizon = 50;
  double safety_distance = 30.0;

  h.onMessage([&map_x, &map_y, &map_s, &map_dx, map_dy, &ref_lane, &ref_vel, &target_vel, &safety_distance, &horizon](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();


        if (event == "telemetry") {
            // j[1] is the data JSON object

            // Main car's localization Data
            double car_x = j[1]["x"];
            double car_y = j[1]["y"];
            double car_s = j[1]["s"];
            double car_d = j[1]["d"];
            double car_yaw = j[1]["yaw"];
            double car_speed = j[1]["speed"];

            // Previous path data given to the Planner
            auto previous_path_x = j[1]["previous_path_x"];
            auto previous_path_y = j[1]["previous_path_y"];
            // Previous path's end s and d values
            double end_path_s = j[1]["end_path_s"];
            double end_path_d = j[1]["end_path_d"];

            // Sensor Fusion Data, a list of all other cars on the same side of the road.
            vector<vector<double>> sensor_fusion = j[1]["sensor_fusion"];

            int prev_size = previous_path_x.size();

            static int counter = 0;



            vector<double> ptsx, ptsy;

            double ref_x = car_x;
            double ref_y = car_y;
            double ref_yaw = deg2rad(car_yaw);
            double ref_s = prev_size > 0 ? end_path_s : car_s;
            bool too_close = false;
            double close_object_speed;
            double close_object_s = numeric_limits<double>::max();


            if(++counter % 50 == 0) {

                double min_s = numeric_limits<double>::max();
                double max_s = numeric_limits<double>::min();
                for(int i = 0; i < sensor_fusion.size(); ++i) {
                    double s = sensor_fusion[i][5];
                    min_s = std::min(min_s, s);
                    max_s = std::max(max_s, s);
                }
                cout << "sensor s: " << min_s << "  " << car_s << "  " << max_s << endl;
                vector<double> lane_speeds(3, target_vel);
                caluclateEffectiveLaneSpeeds(lane_speeds, sensor_fusion, ref_s, prev_size/50.0, safety_distance);
                char buffer[300];
                sprintf(buffer,"lane speeds: %4.2f  |  %4.2f  |  %4.2f", lane_speeds[0], lane_speeds[1], lane_speeds[2]);
                cout << buffer << endl;

                vector<double> lane_costs;
                calculateEfficiencyCost(lane_costs, lane_speeds, ref_lane, target_vel);
                memset(buffer, 300, 0);
                sprintf(buffer,"lane costs:  %4.2f  |  %4.2f  |  %4.2f", lane_costs[0], lane_costs[1], lane_costs[2]);
                cout << buffer << endl;

                vector<double> safety_costs;
                calculateSafetyCost(safety_costs, sensor_fusion, ref_lane, ref_s, prev_size/50.0, safety_distance);
                memset(buffer, 300, 0);
                sprintf(buffer,"safe costs:  %4.2f  |  %4.2f  |  %4.2f", safety_costs[0], safety_costs[1], safety_costs[2]);
                cout << buffer << endl;

                int min_cost_lane = ref_lane;
                double min_cost = numeric_limits<double>::max();
                for(int i = 0; i < 3; ++i) {
                    double c = lane_costs[i] + 1000 * safety_costs[i];
                    if(c < min_cost) {
                        min_cost = c;
                        min_cost_lane = i;
                    }
                }

                if(min_cost_lane != ref_lane) {
                    ref_lane = min_cost_lane;
                    cout << "changing lane to " << ref_lane << endl;
                }
                cout << endl;
            }

            for(int object_id = 0; object_id < sensor_fusion.size(); ++object_id) {
                auto& object = sensor_fusion[object_id];
                double object_d = object[6];
                int object_lane = static_cast<int>((object_d - 2.0) / 4.0 + 0.5);

                if(object_lane == ref_lane) {
                    double object_s = object[5];
                    double object_vx = object[3];
                    double object_vy = object[4];
                    double object_speed = sqrt(object_vx*object_vx + object_vy*object_vy);
                    object_s += prev_size * 0.02 * object_speed;

                    if(object_s > ref_s) {
                        double object_distance = object_s - ref_s;
                        //cout << "object distance = " << object_distance << endl;
                        if(object_distance < safety_distance) {
                            too_close = true;
                            if(object_s < close_object_s) {
                                close_object_s = object_s;
                                close_object_speed = object_speed;
                            }
                        }

                    }
                }
            }




            if(prev_size < 2) {
                double prev_car_x = car_x - cos(car_yaw);
                double prev_car_y = car_y - sin(car_yaw);

                ptsx.push_back(prev_car_x);
                ptsx.push_back(car_x);

                ptsy.push_back(prev_car_y);
                ptsy.push_back(car_y);
            } else {

                ref_x = previous_path_x.back();
                ref_y = previous_path_y.back();

                double ref_x_prev = previous_path_x[prev_size-2];
                double ref_y_prev = previous_path_y[prev_size-2];

                ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

                ptsx.push_back(ref_x_prev);
                ptsx.push_back(ref_x);

                ptsy.push_back(ref_y_prev);
                ptsy.push_back(ref_y);
            }


            double wp_step_size = 30.0;
            for(int i = 1; i <= 3; ++i) {
                auto xy = getXY(ref_s + i*wp_step_size, (2 + 4 * ref_lane), map_s, map_x, map_y);
                ptsx.push_back(xy[0]);
                ptsy.push_back(xy[1]);
            }

            for(int i = 0; i < ptsx.size(); ++i) {
                double shift_x = ptsx[i] - ref_x;
                double shift_y = ptsy[i] - ref_y;

                ptsx[i] = (shift_x * cos(0-ref_yaw) - shift_y * sin(0-ref_yaw));
                ptsy[i] = (shift_x * sin(0-ref_yaw) + shift_y * cos(0-ref_yaw));
            }

            tk::spline wp_spline;
            wp_spline.set_points(ptsx, ptsy);

            vector<double> next_x_vals(prev_size), next_y_vals(prev_size);

            std::copy(previous_path_x.begin(), previous_path_x.end(), next_x_vals.begin());
            std::copy(previous_path_y.begin(), previous_path_y.end(), next_y_vals.begin());

            double target_x = wp_step_size;
            double target_y = wp_spline(target_x);
            double target_dist = sqrt(target_x*target_x + target_y*target_y);

            double x_add_on = 0.0;

            double projected_ref_s = ref_s;
            for(int i = 0; i < horizon - prev_size; ++i) {
                if(too_close) {
                    double projected_object_s = close_object_s + (prev_size + i) * 0.02 * close_object_speed;
                    if(projected_ref_s > projected_object_s - safety_distance)
                        ref_vel = std::max(ref_vel - 0.3, 0.0);
                    else
                        ref_vel = std::min(ref_vel + 0.1, target_vel);
                } else {
                    ref_vel = std::min(ref_vel + 0.224, target_vel);
                }

                double frame_vel = 0.02 * mphToMs(ref_vel);
                projected_ref_s += frame_vel;
                double N = target_dist / frame_vel;
                double x_point = x_add_on + target_x / N;
                double y_point = wp_spline(x_point);

                x_add_on = x_point;

                double x_ref = x_point;
                double y_ref = y_point;

                x_point = (x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw));
                y_point = (x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw));

                x_point += ref_x;
                y_point += ref_y;

                next_x_vals.push_back(x_point);
                next_y_vals.push_back(y_point);
            }

            json msgJson;
            // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
            msgJson["next_x"] = next_x_vals;
            msgJson["next_y"] = next_y_vals;

            auto msg = "42[\"control\","+ msgJson.dump()+"]";

            //this_thread::sleep_for(chrono::milliseconds(1000));
            ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
