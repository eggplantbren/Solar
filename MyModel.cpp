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
    C = -10.0 + 20.0*rng.rand();
    A = 5.0*rng.rand();
    T = 365.0 + 0.5*rng.rand();
    phi = 2.0*M_PI*rng.rand();
    mu = -5.0 + 5.0*rng.rand();
    beta = exp(log(1E-3) + log(1E3)*rng.rand());
    L = exp(log(0.1) + log(1E3)*rng.rand());
}

double MyModel::perturb(DNest4::RNG& rng)
{
    double logH = 0.0;

    int which = rng.rand_int(7);

    if(which == 0)
    {
        C += 20.0*rng.randh();
        DNest4::wrap(C, -10.0, 10.0);
    }
    else if(which == 1)
    {
        A += 5.0*rng.randh();
        DNest4::wrap(A, 0.0, 5.0);
    }
    else if(which == 2)
    {
        T += 10.0*rng.randh();
        DNest4::wrap(T, 365.0, 365.5);
    }
    else if(which == 3)
    {
        phi += 2.0*M_PI*rng.randh();
        DNest4::wrap(phi, 0.0, 2.0*M_PI);
    }
    else if(which == 4)
    {
        mu += 5.0*rng.randh();
        DNest4::wrap(mu, -5.0, 0.0);
    }
    else if(which == 5)
    {
        beta = log(beta);
        beta += log(1E3)*rng.randh();
        DNest4::wrap(beta, log(1E-3), log(1.0));
        beta = exp(beta);
    }
    else
    {
        L = log(L);
        L += log(1E3)*rng.randh();
        DNest4::wrap(L, log(0.1), log(1E2));
        L = exp(L);
    }

    return logH;
}

double MyModel::log_likelihood() const
{
    double logL = 0.0;

    std::vector<double> fracs(y.size());
    std::vector<double> logits(y.size());
    for(size_t i=0; i<y.size(); ++i)
    {
        double top = exp(C + A*sin(2.0*M_PI*t[i]/T + phi));
        fracs[i] = y[i]/top;
        logits[i] = log(fracs[i]/(1.0 - fracs[i]));
    }

    // sigma^2 = a^2*sigma^2 + b^2
    // sigma^2(1 - a^2) = b^2
    double alpha = exp(-1.0/L);
    double sigma = beta/sqrt(1.0 - alpha*alpha);
    logL += -0.5*log(2.0*M_PI*sigma*sigma)
            -0.5*pow(logits[0] - mu, 2)/(sigma*sigma)
            - log(y[0]) - log(1.0 - fracs[0]);  // Jacobian
    for(size_t i=1; i<logits.size(); ++i)
    {
        double expected = mu + alpha*(logits[i-1] - mu);
        logL += -0.5*log(2.0*M_PI*beta*beta)
                -0.5*pow(logits[i] - expected, 2)/(beta*beta)
                - log(y[i]) - log(1.0 - fracs[i]);  // Jacobian
    }

    return logL;
}

void MyModel::print(std::ostream& out) const
{
    out << std::setprecision(12);
    out << C << ' ' << A << ' ' << T << ' ' << phi << ' ';
    out << mu << ' ' << beta << ' ' << L;
}

std::string MyModel::description() const
{
    return std::string("C A T phi mu beta L");
}

