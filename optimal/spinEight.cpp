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


//构造残差函数和雅可比矩阵
class AnalyticCostFunctor : public ceres::SizedCostFunction<1,4> {
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
            const double alpha = parameters[0][3];
            // cout << "theta: " << theta << "ro: " << ro << endl;
            // cout << "ro: " << ro << endl;  
            // 残差矩阵
            

            //场地是一个100*40的有限场地，因此存在约束条件。
            // 1) r + sqrt(2)*r + c < 100, 2) c - r - sqrt(2)*r > 0
            // 3) d + r < 40, 4) d - r > 0
            // 采用惩罚函数法，当不满足约束条件是残差十分大。
            // 判断是曲线还是直线
            if(r + sqrt(2)*r + c < 100 && c - r - sqrt(2)*r > 0 && d + r < 40 && d - r > 0)
            {
              //如果在圆心右侧
              if(data_x > c)
              {
                
                double theta1 = atan2(data_y -d + sqrt(2)*r*sin(alpha), data_x - sqrt(2)*r*cos(alpha) - c);
                //残差矩阵
                double x = (sqrt(2)*r + r*cos(theta1))*cos(alpha) + r*sin(theta1)*sin(alpha) + c;
                double y = -(cos(theta1) + sqrt(2))*r*sin(alpha) + r*sin(theta1)*cos(alpha) + d;
                residuals[0] = pow(data_y - y, 2) + pow(data_x - x, 2);
                if (!jacobians) return true;
                double*jacobian = jacobians[0]; 
                if (!jacobian) return true;
                // 雅可比矩阵
                jacobian[0] = -2*(data_y - y)*(-sin(alpha)*(cos(theta1) + sqrt(2)) + sin(theta1)*cos(alpha))
                            -2*(data_x - x)*(cos(alpha)*(cos(theta1) + sqrt(2)) + sin(theta1)*sin(alpha));
                jacobian[1] = -2*(data_x - x);
                jacobian[2] = -2*(data_y - y);
                jacobian[3] = -2*(data_y - y)*(-sin(alpha)*(-sin(alpha)*(r*cos(theta1) + sqrt(2)*r) + r*sin(theta1)*cos(alpha)))
                                -2*(data_x - x)*(-cos(alpha)*(cos(theta1) + sqrt(2))*r - r*sin(theta1)*sin(alpha));
              
              }
              else
              {
                double theta1 = atan2(data_y -d - sqrt(2)*r*sin(alpha), data_x + sqrt(2)*r*cos(alpha) - c);
                //残差矩阵
                double x = (-sqrt(2)*r + r*cos(theta1))*cos(alpha) + r*sin(theta1)*sin(alpha) + c;
                double y = -(cos(theta1) - sqrt(2))*r*sin(alpha) + r*sin(theta1)*cos(alpha) + d;

                residuals[0] = pow(data_y - y, 2) + pow(data_x - x, 2);
                if (!jacobians) return true;
                double*jacobian = jacobians[0];
                if (!jacobian) return true;
                // 雅可比矩阵
                jacobian[0] = -2*(data_y - y)*(-sin(alpha)*(cos(theta1) - sqrt(2)) + sin(theta1)*cos(alpha))
                            -2*(data_x - x)*(cos(alpha)*(cos(theta1) - sqrt(2)) + sin(theta1)*sin(alpha));
                jacobian[1] = -2*(data_x - x);
                jacobian[2] = -2*(data_y - y);
                jacobian[3] = -2*(data_y - y)*(-sin(alpha)*(-sin(alpha)*(r*cos(theta1) - sqrt(2)*r) + r*sin(theta1)*cos(alpha)))
                                -2*(data_x - x)*(-cos(alpha)*(cos(theta1) - sqrt(2))*r - r*sin(theta1)*sin(alpha));
              }
            }
          else
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



 int main()
 {

    //直线加圆弧轨迹构建
    double speed = 2;
    double r = 10;
    double c = 45;
    double d = 20;
    double alpha = 0;
    double sampleTime = 0.1;
    vector<double> data_x;
    vector<double> data_y;
    double T = 2*r/speed;
    for(double t = 0; t < T; t=t+0.1)
    {
      double x = -sqrt(2)/2*r + sqrt(2)/2*speed*t;
      double y = -sqrt(2)/2*r + sqrt(2)/2*speed*t;
      data_x.push_back(cos(alpha)*x + sin(alpha)*y + c);
      data_y.push_back(-sin(alpha)*x + cos(alpha)*y + d);
    }

    T = 3*pi*r/2/speed;
    for(double t = 0; t < T; t=t+0.1)
    {
      double theta1 = 3*pi/4 - speed*t/r;
      double x = sqrt(2)*r + r*cos(theta1);
      double y = r*sin(theta1);
      data_x.push_back(x*cos(alpha) + y*sin(alpha) + c);
      data_y.push_back(-x*sin(alpha) + y*cos(alpha) + d);
    }

    T = 2*r/speed;
    for(double t = 0; t < T; t=t+0.1)
    {
      double theta1 = 3*pi/4 - speed*t/r;
      double x = sqrt(2)/2*r - sqrt(2)/2*speed*t;
      double y = -sqrt(2)/2*r + sqrt(2)/2*speed*t;
      data_x.push_back(x*cos(alpha) + y*sin(alpha) + c);
      data_y.push_back(-x*sin(alpha) + y*cos(alpha) + d);     
    }

    T = 3*pi*r/2/speed;
    for(double t = 0; t < T; t=t+0.1)
    {
      double theta2 = pi/4 + speed*t/r;
      double x = -sqrt(2)*r + r*cos(theta2);
      double y = r*sin(theta2);
      data_x.push_back(x*cos(alpha) + y*sin(alpha) + c);
      data_y.push_back(-x*sin(alpha) + y*cos(alpha) + d);
    }

    double theta;
    double ro;
    //构造8字形函数
    double a_init = 10;
    double b_init = 40;
    double c_init = 16.5;
    double d_init = 0;
    double param[4] = {a_init,b_init,c_init,d_init};


    Problem problem;
    for (int i = data_x.size()/3; i < data_x.size()/2; i++) {
    // 对每个点加随机值
    // CostFunction *cost_function = new AnalyticCostFunctor(data_x[i]+(float)(rand()%20)/20-1.0, data_y[i]+(float)(rand()%20)/20-1.0);
    CostFunction *cost_function = new AnalyticCostFunctor(data_x[i], data_y[i]);
    problem.AddResidualBlock(cost_function, NULL, param);
    }

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
    std::cout << "d : " << d_init
            << " -> " << param[3] << "\n";        
    cv::Mat figure(1000, 1000, CV_8UC1,cv::Scalar(255) );;

    plot(figure,r,c,d,alpha);
    plot(figure,param[0],param[1],param[2],param[3]);
    cv::imshow("figure",figure);
    cv::waitKey(0);
    //  std::cout << "d : " << d_init
    //             << " -> " << param[3] << "\n";
    return 0;




}