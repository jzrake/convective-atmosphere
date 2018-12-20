#pragma once
#include <cmath>
#include <algorithm>




// ============================================================================
namespace newtonian_hydro {

    // Indexes to primitive quanitites P
    enum {
        RHO = 0,
        V11 = 1,
        V22 = 2,
        V33 = 3,
        PRE = 4,
    };

    // Indexes to conserved quanitites U
    enum {
        DDD = 0,
        S11 = 1,
        S22 = 2,
        S33 = 3,
        NRG = 4,
    };

    using Vars = std::array<double, 5>;
    using Unit = std::array<double, 3>;
    using Position = std::array<double, 2>;

    struct cons_to_prim;
    struct prim_to_cons;
    struct prim_to_flux;
    struct prim_to_eval;
    struct riemann_hlle;
    struct source_terms;
}




// ============================================================================
struct newtonian_hydro::cons_to_prim
{
    inline Vars operator()(Vars U) const
    {
        const double gm1 = gammaLawIndex - 1.0;
        const double pp = U[S11] * U[S11] + U[S22] * U[S22] + U[S33] * U[S33];
        auto P = Vars();

        P[RHO] =  U[DDD];
        P[PRE] = (U[NRG] - 0.5 * pp / U[DDD]) * gm1;
        P[V11] =  U[S11] / U[DDD];
        P[V22] =  U[S22] / U[DDD];
        P[V33] =  U[S33] / U[DDD];

        return P;
    }
    double gammaLawIndex = 5. / 3;
};




// ============================================================================
struct newtonian_hydro::prim_to_cons
{
    inline Vars operator()(Vars P) const
    {
        const double gm1 = gammaLawIndex - 1.0;
        const double vv = P[V11] * P[V11] + P[V22] * P[V22] + P[V33] * P[V33];
        auto U = Vars();

        U[DDD] = P[RHO];
        U[S11] = P[RHO] * P[V11];
        U[S22] = P[RHO] * P[V22];
        U[S33] = P[RHO] * P[V33];
        U[NRG] = P[RHO] * 0.5 * vv + P[PRE] / gm1;

        return U;
    }
    double gammaLawIndex = 5. / 3;
};




// ============================================================================
struct newtonian_hydro::prim_to_flux
{
    inline Vars operator()(Vars P, Unit N) const
    {
        const double vn = P[V11] * N[0] + P[V22] * N[1] + P[V33] * N[2];
        auto U = prim_to_cons()(P);
        auto F = Vars();

        F[DDD] = vn * U[DDD];
        F[S11] = vn * U[S11] + P[PRE] * N[0];
        F[S22] = vn * U[S22] + P[PRE] * N[1];
        F[S33] = vn * U[S33] + P[PRE] * N[2];
        F[NRG] = vn * U[NRG] + P[PRE] * vn;

        return F;
    }
    double gammaLawIndex = 5. / 3;
};




// ============================================================================
struct newtonian_hydro::prim_to_eval
{
    inline Vars operator()(Vars P, Unit N) const
    {
        const double gm0 = gammaLawIndex;
        const double cs = std::sqrt(gm0 * P[PRE] / P[RHO]);
        const double vn = P[V11] * N[0] + P[V22] * N[1] + P[V33] * N[2];
        auto A = Vars();

        A[0] = vn - cs;
        A[1] = vn;
        A[2] = vn;
        A[3] = vn;
        A[4] = vn + cs;

        return A;
    }
    double gammaLawIndex = 5. / 3;
};




// ============================================================================
struct newtonian_hydro::riemann_hlle
{
    riemann_hlle(Unit nhat) : nhat(nhat) {}

    inline Vars operator()(Vars Pl, Vars Pr) const
    {
        auto Ul = p2c(Pl);
        auto Ur = p2c(Pr);
        auto Al = p2a(Pl, nhat);
        auto Ar = p2a(Pr, nhat);
        auto Fl = p2f(Pl, nhat);
        auto Fr = p2f(Pr, nhat);

        const double epl = *std::max_element(Al.begin(), Al.end());
        const double eml = *std::min_element(Al.begin(), Al.end());
        const double epr = *std::max_element(Ar.begin(), Ar.end());
        const double emr = *std::min_element(Ar.begin(), Ar.end());
        const double ap = std::max(0.0, std::max(epl, epr));
        const double am = std::min(0.0, std::min(eml, emr));

        Vars U, F;

        for (int q = 0; q < 5; ++q)
        {
            U[q] = (ap * Ur[q] - am * Ul[q] + (Fl[q] - Fr[q])) / (ap - am);
            F[q] = (ap * Fl[q] - am * Fr[q] - (Ul[q] - Ur[q]) * ap * am) / (ap - am);
        }
        return F;
    }
    Unit nhat;
    prim_to_cons p2c;
    prim_to_eval p2a;
    prim_to_flux p2f;
};




// ============================================================================
struct newtonian_hydro::source_terms
{
    inline Vars operator()(Vars P, Position X) const
    {
        const double r = X[0];
        const double q = X[1];
        const double dg = P[0];
        const double vr = P[1];
        const double vq = P[2];
        const double vp = P[3];
        const double pg = P[4];
        auto S = Vars();


        // Source terms for spherical geometry.
        // --------------------------------------------------------------------
        S[DDD] = 0.0;
        S[S11] = (2 * pg + dg * (vq * vq + vp * vp)) / r;
        S[S22] = (pg * cot(q) + dg * (vp * vp * cot(q) - vr * vq)) / r;
        S[S33] = -dg * vp * (vr + vq * cot(q)) / r;
        S[NRG] = 0.0;


        // Source terms for point mass gravity. GM = 1.0.
        // --------------------------------------------------------------------
        const double g = 1.0 / r / r;
        S[S11] += -dg * g;
        S[NRG] += -dg * g * vr;


        // Source terms for thermal heating
        // --------------------------------------------------------------------
        S[NRG] += std::exp(-r * r) * std::sin(q * 6) * 0.25;


        return S;
    }
    double cot(double x) const
    {
        return std::tan(M_PI_2 - x);
    }
};
