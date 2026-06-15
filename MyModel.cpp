#include "MyModel.h"
#include "DNest4/code/DNest4.h"
#include <fstream>

std::vector<double> MyModel::t;
std::vector<double> MyModel::y;
std::vector<double> MyModel::log_y;

MyModel::MyModel()
{
    if(t.size() == 0)
    {
        std::fstream fin("data.txt", std::ios::in);
        double _t, _y;
        while(fin >> _t && fin >> _y)
        {
            t.push_back(_t);
            y.push_back(_y);
            log_y.push_back(log(_y));
        }
        fin.close();
    }
}

void MyModel::from_prior(DNest4::RNG& rng)
{
    // Lognormal
    A = exp(10.0*rng.randn());
    T = 360.0 + 10.0*rng.rand();
    phi = 2.0*M_PI*rng.rand();
    L = 5.0*rng.rand(); // Attenuation on log scale
}

double MyModel::perturb(DNest4::RNG& rng)
{
    double logH = 0.0;

    int which = rng.rand_int(4);

    if(which == 0)
    {
        A = log(A);
        logH -= -0.5*pow(A/10.0, 2);
        A += 10.0*rng.randh();
        logH += -0.5*pow(A/10.0, 2);
        A = exp(A);
    }
    else if(which == 1)
    {
        T += 10.0*rng.randh();
        DNest4::wrap(T, 360.0, 370.0);
    }
    else if(which == 2)
    {
        phi += 2.0*M_PI*rng.randh();
        DNest4::wrap(phi, 0.0, 2.0*M_PI);
    }
    else
    {
        L += 5.0*rng.randh();
        DNest4::wrap(L, 0.0, 5.0);
    }

    return logH;
}

double MyModel::log_likelihood() const
{
    double logL = 0.0;

    for(size_t i=0; i<y.size(); ++i)
    {
        double top = A*sin(2.0*M_PI*t[i]/T + phi);

        if(log_y[i] > top)
            logL += -1E300;

        logL += -log(L) + (log_y[i] - top)/L;
    }

    return logL;
}

void MyModel::print(std::ostream& out) const
{
    out << std::setprecision(12);
    out << A << ' ' << T << ' ' << phi << ' ' << L;
}

std::string MyModel::description() const
{
    return std::string("A T phi L");
}

