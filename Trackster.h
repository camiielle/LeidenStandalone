#ifndef TRACKSTER_H
#define TRACKSTER_H

#include <TVector3.h>
#include <iostream>

class Barycenter {
  public:
    Barycenter() = default;
    Barycenter(float eta, float phi, float x, float y, float z) : eta_(eta), phi_(phi), x_(x), y_(y), z_(z) {};
    float eta() const {return eta_;};
    float phi() const {return phi_;};
    float x() const {return x_;};
    float y() const {return y_;};
    float z() const {return z_;};

    void setBarycenter(float eta, float phi, float x, float y, float z){
      eta_ = eta;
      phi_ = phi;
      x_ = x;
      y_ = y;
      z_ = z;
    }



  private:
    float eta_;
    float phi_;
    float x_;
    float y_;
    float z_;
};

class Trackster {
public:

    Trackster() : barycenter_(0,0,0,0,0), raw_energy_(0) {}
    
    void setBarycenter(float eta, float phi, float x, float y, float z){
      barycenter_.setBarycenter(eta,phi,x,y,z);
    }

    const Barycenter barycenter() const {
        return barycenter_;
    }

    void setRawEnergy(float raw_energy) {
      raw_energy_ = raw_energy;
    }

    const float raw_energy() const {
        return raw_energy_;
    }

    void Print() const {
        std::cout << "Trackster: "
                  << "barycenter (x, y, z) = (" 
                  << barycenter().x() << ", " 
                  << barycenter().y() << ", " 
                  << barycenter().z() << "), "
                  << "barycenter_eta = " << barycenter().eta() << ", "
                  << "barycenter_phi = " << barycenter().phi() << ", "
                  << "raw_energy = " << raw_energy() << std::endl;
    }
private:
    Barycenter barycenter_;
    float raw_energy_;
};

#endif // TRACKSTER_H

