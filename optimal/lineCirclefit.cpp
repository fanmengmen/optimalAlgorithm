#include<iostream>
#include<ceres/ceres.h>
#include<math.h>
#include<vector>
#include "plot.hpp"
using namespace std;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
#define pi 3.141593
#include <random>


class AnalyticCostFunctor : public ceres::SizedCostFunction<1,3> {
   public:
     AnalyticCostFunctor(const double x, const double y) : data_x(x), data_y(y) {}
     virtual ~AnalyticCostFunctor() {}
     virtual bool Evaluate(double const* const* parameters,
                           double* residuals,
                           double** jacobians) const {
            // const double speed = parameters[0][0];
            const double r = parameters[0][0];
            const double c = parameters[0][1];
            const double d = parameters[0][2];
            // cout << "theta: " << theta << "ro: " << ro << endl;
            // cout << "ro: " << ro << endl;  
            
            /***************场地是一个100*40的有限场地，因此存在约束条件。
                1) r + sqrt(2)*r + c < 100, 2) c - r - sqrt(2)*r > 0
                3) d + r < 40, 4) d - r > 0
                采用惩罚函数法，当不满足约束条件是残差十分大。
                判断是曲线还是直线(未包含)***************************/
            if(r + sqrt(2)*r + c < 100 && c - r - sqrt(2)*r > 0 && d + r < 40 && d - r > 0)
            {
              //如果在圆心右侧
              if(data_x > c)
              {
                double theta1 = atan2(data_y - d, data_x - sqrt(2)*r - c);
                //残差矩阵
                residuals[0] = pow(data_y - d - r*sin(theta1) ,2) + pow(data_x - r*cos(theta1) - sqrt(2)*r - c,2);
                if (!jacobians) return true;
                double*jacobian = jacobians[0];
                if (!jacobian) return true;
                // 雅可比矩阵
                jacobian[0] = -sin(theta1)*(data_y - d - r*sin(theta1))*2 - (cos(theta1) + sqrt(2))*(data_x - r*cos(theta1) - sqrt(2)*r - c)*2;
                jacobian[1] = -(data_x - r*cos(theta1) - sqrt(2)*r - c)*2;
                jacobian[2] = -(data_y - d - r*sin(theta1))*2;
              
              }
              else//圆心在左侧
              {
                double theta2 = atan2(data_y - d, data_x + sqrt(2)*r - c);
                residuals[0] = pow(data_y - d - r*sin(theta2) ,2) + pow(data_x - r*cos(theta2) + sqrt(2)*r - c,2);
                if (!jacobians) return true;
                double*jacobian = jacobians[0];
                if (!jacobian) return true;
                // 雅可比矩阵
                jacobian[0] = -sin(theta2)*(data_y - d - r*sin(theta2))*2 - cos(theta2) - sqrt(2);
                jacobian[1] = -(data_x - r*cos(theta2) + sqrt(2)*r - c)*2;
                jacobian[2] = -(data_y - d - r*sin(theta2))*2;
              }
            }
            else//不满足约束条件
            {
              residuals[0] = 1000;
              residuals[1] = 1000;
              residuals[2] = 1000;
             
            }
             
            return true;
     }

   private:
     const double data_x;
     const double data_y;
 };





int main(int argc, char** argv) {
    // google::InitGoogleLogging(argv[0]);

   
  

    srand(time(NULL));
    //直线加圆弧轨迹构建
    double speed = 2;
    double r = 10;
    double c = 45;
    double d = 20;
    double sampleTime = 0.1;
    vector<double> data_x;
    vector<double> data_y;
    double T = 2*r/speed;
    for(double t = 0; t < T; t=t+0.1)
    {
      data_x.push_back(-sqrt(2)/2*r + sqrt(2)/2*speed*t + c);
      data_y.push_back(-sqrt(2)/2*r + sqrt(2)/2*speed*t + d);
    }
    T = 3*pi*r/2/speed;
    for(double t = 0; t < T; t=t+0.1)
    {
      double theta1 = 3*pi/4 - speed*t/r;
      data_x.push_back(sqrt(2)*r + r*cos(theta1) + c);
      data_y.push_back(r*sin(theta1) + d);
    }
    T = 2*r/speed;
     for(double t = 0; t < T; t=t+0.1)
    {
      double theta1 = 3*pi/4 - speed*t/r;
      data_x.push_back(sqrt(2)/2*r - sqrt(2)/2*speed*t + c);
      data_y.push_back(-sqrt(2)/2*r + sqrt(2)/2*speed*t + d);      
    }
    T = 3*pi*r/2/speed;
    for(double t = 0; t < T; t=t+0.1)
    {
      double theta2 = pi/4 + speed*t/r;
      data_x.push_back(-sqrt(2)*r + r*cos(theta2) + c);
      data_y.push_back(r*sin(theta2) + d);
    }


    // 设置初始化参数
    double a_init = 15;
    double b_init = 40;
    double c_init = 16.5;
    double param[3] = {a_init,b_init,c_init};
    // cout << data[1] << "   " << data[2] << endl;

    //检测初值是否满足约束条件，且有意义
    // cout << "go_____: " << data[1]-b*ro*sin(theta) + data[0]-a*ro*cos(theta) << endl;


    // Build the problem.
    Problem problem;
    // 加随机值
    for (int i = data_x.size()/3; i < data_x.size()/2; i++) {
        CostFunction *cost_function = new AnalyticCostFunctor(data_x[i]+(float)(rand()%20)/20-1.0,
             data_y[i]+(float)(rand()%20)/20-1.0);
        problem.AddResidualBlock(cost_function, NULL, param);
    }
    // 不加随机值
    // for (int i = data_x.size()/3; i < data_x.size()*2/5; i++) {
    //     CostFunction *cost_function = new AnalyticCostFunctor(data_x[i],
    //          data_y[i]);
    //     problem.AddResidualBlock(cost_function, NULL, param);
    // }
    
    // Run the solver!
    Solver::Options options;
    options.gradient_tolerance = 1e-32;
    options.function_tolerance = 1e-16;
    options.linear_solver_type = ceres::DENSE_QR;
    options.max_num_iterations = 100;
    options.minimizer_progress_to_stdout = true;
    Solver::Summary summary;
    Solve(options, &problem, &summary);
    std::cout << summary.BriefReport() << "\n";
    std::cout << "a : " << a_init
                << " -> " << param[0] << "\n";
    std::cout << "b : " << b_init
                << " -> " << param[1] << "\n";
    std::cout << "c : " << c_init
                << " -> " << param[2] << "\n";

    cv::Mat figure(1000, 1000, CV_8UC1,cv::Scalar(255) );
    plot(figure,r,c,d);
    plot(figure,param[0],param[1],param[2]);
    cv::imshow("figure",figure);
    cv::waitKey(0);
    //  std::cout << "d : " << d_init
    //             << " -> " << param[3] << "\n";
    return 0;




}


