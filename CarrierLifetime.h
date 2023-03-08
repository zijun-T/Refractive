#pragma once


class Hello
{
    public:
        double lifetime(double rho)
    {
        double t = 0;
        double r, p1, p2, p3, p4, p5;
        r = rho * 1e6;
        //cout << r <<endl;
        
        if (Element == "Si")
        {
            p1 = -0.001098;
            p2 = 0.05413;
            p3 = -1.03;
            p4 = 8.811;
            p5 = -31.2;
            
            if (rho < 1e12)
            {
                r = 1e12 * 1e6;
            }
            else if(rho > 1e20)
            {
                r = 1e20 * 1e6;
            }
        }
        else
        {
            cout << Element + "lifetime loss!!!" << endl;
            exit(0);
        }

        r = log10(r);
        t = p1*r*r*r*r + p2*r*r*r + p3*r*r + p4*r +p5;
        t = pow(10,t);

        //cout << t <<endl;
        //exit(0);
        return t;
    }
};




















