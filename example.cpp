/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2014 Vlad Gheorghiu (vgheorgh@gmail.com)
 *
 * This file is part of Quantum++.
 *
 * Quantum++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Quantum++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Quantum++.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "qpp.h"

// #include "MATLAB/matlab.h" // support for MATLAB

using namespace qpp;

int main()
{
	//	/*   comment this line with "//" to uncomment the whole example

	// Qudit Teleportation
	{
		std::size_t D = 3; // size of the system
		std::cout << std::endl << "**** Qudit Teleportation, D = " << D
				<< " ****" << std::endl;
		ket mes_AB = ket::Zero(D * D); // maximally entangled state resource
		for (std::size_t i = 0; i < D; i++)
			mes_AB += mket( { i, i }, { D, D });
		mes_AB /= std::sqrt((double) D);
		cmat Bell_aA = adjoint( // circuit that measures in the qudit Bell basis
				gt.CTRL(gt.Xd(D), { 0 }, { 1 }, 2, D)
						* kron(gt.Fd(D), gt.Id(D)));
		ket psi_a = randket(D); // random state as input on a
		std::cout << ">> Initial state:" << std::endl;
		displn(psi_a);
		ket input_aAB = kron(psi_a, mes_AB); // joint input state aAB
		// output before measurement
		ket output_aAB = apply(input_aAB, Bell_aA, { 0, 1 }, 3, D);
		auto measured_aA = measure(ptrace2(prj(output_aAB), { D * D, D }),
				gt.Id(D * D)); // measure on aA
		std::discrete_distribution<std::size_t> dd(measured_aA.first.begin(),
				measured_aA.first.end());
		std::cout << ">> Measurement probabilities: ";
		displn(measured_aA.first, ", ");
		std::size_t m = dd(rdevs._rng); // sample
		auto midx = n2multiidx(m, { D, D });
		std::cout << ">> Measurement result: ";
		displn(midx, " ");
		// conditional result on B before correction
		ket output_m_aAB = apply(output_aAB, prj(mket(midx, { D, D })),
				{ 0, 1 }, 3, D) / std::sqrt(measured_aA.first[m]);
		cmat correction_B = powm(gt.Zd(D), midx[0])
				* powm(adjoint(gt.Xd(D)), midx[1]); // correction operator
				// apply correction on B
		output_aAB = apply(output_m_aAB, correction_B, { 2 }, 3, D);
		cmat rho_B = ptrace1(prj(output_aAB), { D * D, D });
		std::cout << ">> Bob's density operator: " << std::endl;
		displn(rho_B);
		std::cout << ">> Norm difference: " << norm<2>(rho_B - prj(psi_a))
				<< std::endl; // verification
	}

	// Qudit Dense Coding
	{
		std::size_t D = 3; // size of the system
		std::cout << std::endl << "**** Qudit Dense Coding, D = " << D
				<< " ****" << std::endl;
		ket mes_AB = ket::Zero(D * D); // maximally entangled state resource
		for (std::size_t i = 0; i < D; i++)
			mes_AB += mket( { i, i }, { D, D });
		mes_AB /= std::sqrt((double) D);
		cmat Bell_AB = adjoint( // circuit that measures in the qudit Bell basis
				gt.CTRL(gt.Xd(D), { 0 }, { 1 }, 2, D)
						* kron(gt.Fd(D), gt.Id(D)));
		// equal probabilities of choosing a message
		std::uniform_int_distribution<std::size_t> uid(0, D * D - 1);
		std::size_t m_A = uid(rdevs._rng); // sample, obtain the message index
		auto midx = n2multiidx(m_A, { D, D });
		std::cout << ">> Alice sent: ";
		displn(midx, " ");
		// Alice's operation
		cmat U_A = powm(gt.Zd(D), midx[0]) * powm(adjoint(gt.Xd(D)), midx[1]);
		// Alice encodes the message
		ket psi_AB = apply(mes_AB, U_A, { 0 }, 2, D);
		// Bob measures the joint system in the qudit Bell basis
		psi_AB = apply(psi_AB, Bell_AB, { 0, 1 }, 2, D);
		auto measured = measure(psi_AB, gt.Id(D * D));
		std::cout << ">> Bob measurement probabilities: ";
		displn(measured.first, ", ");
		// Bob samples according to the measurement probabilities
		std::discrete_distribution<std::size_t> dd(measured.first.begin(),
				measured.first.end());
		std::size_t m_B = dd(rdevs._rng);
		std::cout << ">> Bob received: ";
		displn(n2multiidx(m_B, { D, D }), " ");
	}

	// Entanglement
	{
		std::cout << std::endl << "**** Entanglement ****" << std::endl;

		cmat rho = 0.2 * st.pb00 + 0.8 * st.pb11;
		std::cout << ">> rho: " << std::endl;
		displn(rho);
		std::cout << ">> Concurrence of rho: " << concurrence(rho) << std::endl;
		std::cout << ">> Negativity of rho: " << negativity(rho, { 2, 2 })
				<< std::endl;
		std::cout << ">> Logarithimc negativity of rho: "
				<< lognegativity(rho, { 2, 2 }) << std::endl;
		ket psi = 0.6 * mket( { 0, 0 }) + 0.8 * mket( { 1, 1 });
		// apply some local random unitaries
		psi = kron(randU(2), randU(2)) * psi;
		std::cout << ">> psi: " << std::endl;
		displn(psi);
		std::cout << ">> Entanglement of psi: " << entanglement(psi, { 2, 2 })
				<< std::endl;
		std::cout << ">> Concurrence of psi: " << concurrence(prj(psi))
				<< std::endl;
		std::cout << ">> G-Concurrence of psi: " << gconcurrence(psi)
				<< std::endl;
		std::cout << ">> Schmidt coefficients of psi: " << std::endl;
		displn(schmidtcoeff(psi, { 2, 2 }));
	}

	// */
}
