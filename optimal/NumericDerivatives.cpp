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

// class AnalyticCostFunctor : public ceres::SizedCostFunction<1,4> {
//    public:
//      AnalyticCostFunctor(const double x, const double y) : data_x(x), data_y(y) {}
//      virtual ~AnalyticCostFunctor() {}
//      virtual bool Evaluate(double const* const* parameters,
//                            double* residuals,
//                            double** jacobians) const {
//             const double a = parameters[0][0];
//             const double b = parameters[0][1];
//             const double c = parameters[0][2];
//             const double d = parameters[0][3];
//             const double theta = atan2((data_y-d)*a, (data_x-c)*b);
            
//             double ro;
//             if((theta >= -pi/4 && theta <= pi/4)||theta >= 3*pi/4||theta <= -3*pi/4)
//               ro = sqrt(cos(theta*2));
//             else
//               {
//                 ro = 1;
//               }

//             // cout << "theta: " << theta << "ro: " << ro << endl;
//             // cout << "ro: " << ro << endl;  
//             // 残差矩阵
//             residuals[0] = pow(data_y-b*ro*sin(theta)-d,2) + pow(data_x-a*ro*cos(theta)-c,2);
//             // residuals[0] = pow(data_x - a - b*cos(d),2) + pow(data_y - c - b*sin(d),2);
//             if (!jacobians) return true;
//             double*jacobian = jacobians[0];
//             if (!jacobian) return true;
//             // 雅可比矩阵
//             jacobian[1] = -ro*sin(theta)*(data_y - b*ro*sin(theta)-d)*2;
//             jacobian[0] = -ro*cos(theta)*(data_x - a*ro*cos(theta)-c)*2;
//             jacobian[2] = -2*(data_x - a*ro*cos(theta)-c);
//             jacobian[3] = -2*(data_y - b*ro*sin(theta)-d);

//             // jacobian[0] = -2*(data_x - a - b*cos(d));
//             // jacobian[1] = -cos(d)*(data_x - a - b*cos(d)) -sin(d)*(data_y - c - b*sin(d));
//             // jacobian[2] = -2*(data_y - c - b*sin(d));
//             // jacobian[3] = b*sin(d)*(data_x - a - b*cos(d)) - b*cos(d)*(data_y - c - b*sin(d));

//             return true;
//      }

//    private:
//      const double data_x;
//      const double data_y;
//  };



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
              else
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

// class Rat43Analytic: public SizedCostFunction<1,2>{
//     private:
//         double data_x, data_y;
//         int n;
//         int param_a, param_b;


//     public:
//         Rat43Analytic(const double x, const double y):data_x(x), data_y(y){}
//         virtual ~Rat43Analytic(){}
//         virtual bool Evalute(double const* const* parameters,
//                            double* residual,
//                            double** jacobians)
//         {
//             const double a = parameters[0][0];
//             const double b = parameters[0][1];
//             double theta = atan2(data_y*a, data_x*b);
//             double ro = sqrt(cos(theta*2));
//             residual[0] = pow(data_y-b*ro*sin(theta),2) + pow(data_x-a*ro*cos(theta),2);
//             if (!jacobians) return true;
//             double* jacobian = jacobians[0];
//             if (!jacobian) return true;
//             jacobian[0] = data_x;
//        jacobian[1] = 1.0;
// //        return true;
//             return true;
//         }

        




// };
    ;
//example指数函数损失函数
struct ExponentialResidual {
  ExponentialResidual(double x, double y)
      : x_(x), y_(y) {}

  template <typename T>
  bool operator()(const T* const m, const T* const c, T* residual) const {
    residual[0] = T(y_) - exp(m[0] * T(x_) + c[0]);
    return true;
  }

 private:
  // Observations for a sample.
  const double x_;
  const double y_;
};



int main(int argc, char** argv) {
    // google::InitGoogleLogging(argv[0]);

    //example 指数函数拟合
  
    // int kNumObservations = 108;//有50对点待拟合
    
    // // double data[2*kNumObservations];
    // double a = 40;
    // double b = 40;
    // double c = 40;
    // double d = 30;

    // 构建指数函数，设置参数a和b，将x和y存在data中x[]
    // for(int i=0; i < kNumObservations; i++)
    // {
    //     data[2*i] = 0.1*i;
    //     data[2*i+1] = exp(a*data[2*i] + b);
    //     cout << data[2*i+1] << endl; 
    // }
    // // The variable to solve for with its initial value.
    // double m = 0.0;
    // double c = 10;

    // // Build the problem.
    // Problem problem;
    // for (int i = 0; i < kNumObservations; ++i) {
      // 自动求微分
    // Rat43Analytic* cost_function =
    //     new AutoDiffRat43Analytic<ExponentialResidual, 1, 1, 1>(
    //         new ExponentialResidual(data[2 * i], data[2 * i + 1]));
    // problem.AddResidualBlock(cost_function, NULL, &m, &c);
    // }
    

    // // Run the solver!
    // Solver::Options options;
    // options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY;
    // options.minimizer_progress_to_stdout = true;
    // options.max_linear_solver_iterations = 100;
    // Solver::Summary summary;
    // Solve(options, &problem, &summary);
    // //Show the result and progress!
    // std::cout << summary.BriefReport() << "\n";
    // std::cout << "x : " << 0
    //             << " -> " << m << "\n";
    // std::cout << "x : " << 0
    //             << " -> " << c << "\n";
  
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

    double theta;
    double ro;
    //构造8字形函数
    // for(int i = 0; i < kNumObservations; i++)
    // {
    //     theta = -pi/4 + 0.01*i + 0.01;
    //     ro = sqrt(cos(2*theta));
    //     data[2*i] = a*ro*cos(theta) + c;
    //     data[2*i+1] = b*ro*sin(theta) + d;
    // }

    // The variable to solve for with its initial value.
    // double a_init = 40;
    // double b_init = 35;
    // double c_init = 50;
    // double d_init = 20;
    // double param[4] = {a_init,b_init,c_init,d_init};

    double a_init = 10;
    double b_init = 40;
    double c_init = 16.5;
    double param[3] = {a_init,b_init,c_init};
    // cout << data[1] << "   " << data[2] << endl;
    // cout << "go_____: " << data[1]-b*ro*sin(theta) + data[0]-a*ro*cos(theta) << endl;
    // Build the problem.
    Problem problem;
    for (int i = data_x.size()/3; i < data_x.size()*2/5; i++) {
    CostFunction *cost_function = new AnalyticCostFunctor(data_x[i]+(float)(rand()%20)/20-1.0, data_y[i]+(float)(rand()%20)/20-1.0);
    problem.AddResidualBlock(cost_function, NULL, param);
    }


  // Set up the only cost function (also known as residual). This uses
  // auto-differentiation to obtain the derivative (jacobian).
  
//   // Run the solver!
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


// #include "ceres/ceres.h"
// #include <opencv2/core/core.hpp>
// using namespace std;
// using ceres::CostFunction;
// using ceres::Problem;
// using ceres::Solver;
// class AnalyticCostFunctor : public ceres::SizedCostFunction<1,2> {
//    public:
//      AnalyticCostFunctor(const double x, const double y) : x_(x), y_(y) {}
//      virtual ~AnalyticCostFunctor() {}
//      virtual bool Evaluate(double const* const* parameters,
//                            double* residuals,
//                            double** jacobians) const {
//        const double a = parameters[0][0];
//        const double b = parameters[0][1];

//        residuals[0] = (a*x_+b)-y_;

//        if (!jacobians) return true;
//        double* jacobian = jacobians[0];
//        if (!jacobian) return true;

//        jacobian[0] = x_;
//        jacobian[1] = 1.0;
//        return true;
//      }

//    private:
//      const double x_;
//      const double y_;
//  };


// int main(int argc, char** argv)
// {

//   double a=1.0, b=2.0;         // 真实参数值
//   int N=100;                          // 数据点
//   double w_sigma=0.1;                 // 噪声Sigma值
//   cv::RNG rng;                        // OpenCV随机数产生器
//   double ab[2] = {0,0};            // abc参数的估计值

//   vector<double> x_data, y_data;      // 数据
//   cout<<"generating data: "<<endl;
//   for ( int i=0; i<N; i++ )
//   {
//       double x = i/100.0;
//       x_data.push_back ( x );
//       y_data.push_back (a*x + b + rng.gaussian ( w_sigma ));
//       cout<<x_data[i]<<" "<<y_data[i]<<endl;
//   }
//   // Build the problem.
//   Problem problem;
//   for ( int i=0; i<N; i++ )
//  {
//       CostFunction* cost_function =new AnalyticCostFunctor( x_data[i], y_data[i] );
//       problem.AddResidualBlock(cost_function, NULL, ab);
//   }

//   // Run the solver!
//   Solver::Options options;
//   options.minimizer_progress_to_stdout = true;
//   Solver::Summary summary;
//   Solve(options, &problem, &summary);
//   std::cout << summary.BriefReport() << "\n";
//   std::cout << "output a: " << ab[0]
//             << "output b: " << ab[1] << "\n";
//   return 0;
// }
