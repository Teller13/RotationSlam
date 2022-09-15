#pragma once
#include <memory>  // ???定义了C++标准中的智能指针、内存分配器（allocator）、等函数
#include <Eigen/Dense>
#include <sophus/se3.hpp>
#include <visnav/common_types.h>

namespace visnav{

template <typename Scalar> class AbstractCamera;

template <typename Scalar>
class PinholeCamera : public AbstractCamera<Scalar>
{
  public:
    static constexpr size_t N = 8;
    typedef Eigen::Matrix<Scalar, 2, 1> Vec2;
    typedef Eigen::Matrix<Scalar, 3, 1> Vec3;
    typedef Eigen::Matrix<Scalar, N, 1> VecN;

    PinholeCamera() = default;
    PinholeCamera(const VecN& p) : param(p) {}

    static PinholeCamera<Scalar> getTestProjections() {
      VecN vec1;
      vec1 << 0.5 * 805, 0.5 * 800, 505, 509, 0, 0, 0, 0;
      PinholeCamera<Scalar> res(vec1);

      return res;
    }

    Scalar* data() { return param.data(); }
    const Scalar* data() const { return param.data(); }
    const VecN& getParam() const { return param; }
    static std::string getName() { return "pinhole"; }
    std::string name() const { return getName(); }

    virtual Vec2 project(const Vec3& p) const {
      const Scalar& fx = param[0];
      const Scalar& fy = param[1];
      const Scalar& cx = param[2];
      const Scalar& cy = param[3];

      const Scalar& x = p[0];
      const Scalar& y = p[1];
      const Scalar& z = p[2];

      Vec2 res;

      // TODO SHEET 2: implement camera model
      res(0) = fx * x / z + cx;
      res(1) = fy * y / z + cy;
      return res;
    }
    virtual Vec3 unproject(const Vec2& p) const {
      const Scalar& fx = param[0];
      const Scalar& fy = param[1];
      const Scalar& cx = param[2];
      const Scalar& cy = param[3];

      Vec3 res;

      // TODO SHEET 2: implement camera model
      Vec3 vec;
      vec(0) = (p[0] - cx) / fx;
      vec(1) = (p[0] - cx) / fx;
      vec(2) = Scalar(1);

      res = vec / vec.norm(); //归一化坐标

      return res;
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  private:
    VecN param = VecN::Zero();
};

template <typename Scalar>
class DoubleSphereCamera : public AbstractCamera<Scalar>
{
  public:
    static constexpr size_t N = 8;
    typedef Eigen::Matrix<Scalar, 2, 1> Vec2;
    typedef Eigen::Matrix<Scalar, 3, 1> Vec3;
    typedef Eigen::Matrix<Scalar, N, 1> VecN;

    DoubleSphereCamera() = default;
    DoubleSphereCamera(const VecN& p) : param(p) {}

    static DoubleSphereCamera<Scalar> getTestProjections() {
      VecN vec1;
      vec1 << 0.5 * 805, 0.5 * 800, 505, 509, 0.5 * -0.150694, 0.5 * 1.48785, 0,
          0;
      DoubleSphereCamera<Scalar> res(vec1);

      return res;
    }

    Scalar* data() { return param.data(); }
    const Scalar* data() const { return param.data(); }
    const VecN& getParam() const { return param; }
    static std::string getName() { return "ds"; }
    std::string name() const { return getName(); }


    virtual Vec2 project(const Vec3& p) const {
      const Scalar& fx = param[0];
      const Scalar& fy = param[1];
      const Scalar& cx = param[2];
      const Scalar& cy = param[3];
      const Scalar& xi = param[4];
      const Scalar& alpha = param[5];

      const Scalar& x = p[0];
      const Scalar& y = p[1];
      const Scalar& z = p[2];

      Vec2 res;

      // TODO SHEET 2: implement camera model
      Scalar d1 = p.norm();
      Scalar d2 = sqrt(x * x + y * y + (xi * d1 + z) * (xi * d1 + z));
      res[0] = fx * (x / (alpha * d2 + (Scalar(1) - alpha) * (xi * d1 + z))) + cx;
      res[1] = fy * (y / (alpha * d2 + (Scalar(1) - alpha) * (xi * d1 + z))) + cy;
      return res;
    }
    virtual Vec3 unproject(const Vec2& p) const {
      const Scalar& fx = param[0];
      const Scalar& fy = param[1];
      const Scalar& cx = param[2];
      const Scalar& cy = param[3];
      const Scalar& xi = param[4];
      const Scalar& alpha = param[5];

      Vec3 res;

      // TODO SHEET 2: implement camera model
      Scalar mx = (p[0] - cx) / fx;
      Scalar my = (p[1] - cy) / fy;
      Scalar r2 = mx * mx + my * my;
      Scalar mz =
          ((Scalar)1 - alpha * alpha * r2) /
          (alpha * sqrt((Scalar)1 - ((Scalar)2 * alpha - (Scalar)1) * r2) +
           ((Scalar)1 - alpha));
      Scalar c =
          (mz * xi + sqrt(mz * mz + ((Scalar)1 - xi * xi) * r2)) / (mz * mz + r2);
      res[0] = c * mx;
      res[1] = c * my;
      res[2] = c * mz - xi;

      return res;
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  private:
    VecN param = VecN::Zero();
};

template <typename  Scalar> class AbstractCamera {
  public:
    static constexpr size_t N = 8;  // ***(1) constexpr用法
    typedef Eigen::Matrix<Scalar, 2, 1> Vec2;
    typedef Eigen::Matrix<Scalar, 3, 1> Vec3;
    typedef Eigen::Matrix<Scalar, N, 1> VecN;

    virtual ~AbstractCamera() = default; //***（2）虚函数 (3) 虚析构函数，default用法
    virtual Scalar* data() = 0; // ***(4) 纯虚函数
    virtual const Scalar* data() const = 0;
    virtual Vec2 project(const Vec3& p) const = 0;
    virtual Vec3 unproject(const Vec2& p) const = 0;

    virtual std::string name() const = 0;
    virtual const VecN& getParam() const = 0;

    inline int width() const { return width_; }
    inline int& width() { return width_; }
    inline int height() const { return height_; }
    inline int& height() { return height_; }

    static std::shared_ptr<AbstractCamera> from_data(const std::string& name,
                                                     const Scalar* sIntr
                                                     )
    {
      if (name == DoubleSphereCamera<Scalar>::getName())
      {
        Eigen::Map<Eigen::Matrix<Scalar, 8, 1> const> intr(sIntr);
        return std::shared_ptr<AbstractCamera>(
            new DoubleSphereCamera<Scalar>(intr));
      }
      else if (name == PinholeCamera<Scalar>::getName())
      {
        Eigen::Map<Eigen::Matrix<Scalar, 8, 1> const> intr(sIntr);
        return std::shared_ptr<AbstractCamera>(
            new PinholeCamera<Scalar>(intr));
      }
      else
      {
        std::cerr << "Camera model " << name << " is not implemented."
                  << std::endl;
        std::abort();
      }
    }
    // Loading from double sphere initialization
    static std::shared_ptr<AbstractCamera> initialize(const std::string& name,
                                                      const Scalar* sIntr)
    {
      Eigen::Matrix<Scalar, 8, 1> init_intr;
      if (name == DoubleSphereCamera<Scalar>::getName())
      {
        Eigen::Map<Eigen::Matrix<Scalar, 8, 1> const> intr(sIntr);
        init_intr = intr;
        return std::shared_ptr<AbstractCamera>(
            new DoubleSphereCamera<Scalar>(init_intr));
      }
      else if (name == PinholeCamera<Scalar>::getName())
      {
        Eigen::Map<Eigen::Matrix<Scalar, 8, 1> const> intr(sIntr);
        init_intr = intr;
        init_intr.template tail<4>().setZero();
        return std::shared_ptr<AbstractCamera>(
            new PinholeCamera<Scalar>(init_intr));
      }
      else
      {
        std::cerr << "Camera model " << name << " is not implemented."
                  << std::endl;
        std::abort();
      }
    }

  private:
    // image dimensions
    int width_ = 0;
    int height_ = 0;

};

}

// ***(1) constexpr用法
// int get_five() {return 5;}
// int some_value[get_five() + 7];
// 创建包含12个整数的数组. C++03中非法，因为get_five() + 7不是常量表达式
// Create an array of 12 integers. Valid C++11 合法


// ***(2) 虚函数
// 虚函数最关键的特点是“动态联编”，它可以在运行时判断指针指向的对象，并自动调用相应的函数
// class Parent{
// public:
//    char data[20];
//    void Function1(){这是函数1};
//    virtual void Function2(){这是函数2};
//    }parent;

// ***(3) 虚析构函数
// （1）如果父类的析构函数不加virtual关键字
// 当父类的析构函数不声明成虚析构函数的时候，当子类继承父类，父类的指针指向子类时，
// delete掉父类的指针，只调动父类的析构函数，而不调动子类的析构函数。
// （2）如果父类的析构函数加virtual关键字
// 当父类的析构函数声明成虚析构函数的时候，当子类继承父类，父类的指针指向子类时，
// delete掉父类的指针，先调动子类的析构函数，再调动父类的析构函数。

// ***(4) 纯虚函数
// 纯虚函数是在基类中声明的虚函数，它在基类中没有定义，但要求任何派生类(继承类)都要定义自己的实现方法。
// 在基类中实现纯虚函数的方法是在函数原型后加=0。
// 定义纯虚函数的目的在于，使派生类仅仅只是继承函数的接口。

// ***(5) inline关键字
// 这样可以解决一些频繁调用的函数大量消耗栈空间（栈内存）的问题。
// 关键字inline必须与函数定义放在一起才能使函数成为内联函数，仅仅将inline放在函数声明前面不起任何作用。
// inline是一种“用于实现”的关键字，而不是一种“用于声明”的关键字

// ***(6) shared_ptr智能指针
//在传统 C++ 里我们只好使用 new 和 delete 去 『记得』对资源进行释放。
// 而 C++11 引入智能指针的概念，使用引用计数的想法，让程序员不再需要关心手动释放内存。
//这些智能指针就包括 std::shared_ptr/std::unique_ptr/std::weak_ptr，
// 使用它们需要包含头文件 <memory>

//引用计数不是垃圾回收，引用计数能够尽快收回不再被使用的对象，
// 同时在回收的过程中也不会造成长时间的等待， 更能够清晰明确的表明资源的生命周期。

//std::shared_ptr 是一种智能指针，它能够记录多少个 shared_ptr 共同指向一个对象，
// 从而消除显式的调用 delete，当引用计数变为零的时候就会将对象自动删除。