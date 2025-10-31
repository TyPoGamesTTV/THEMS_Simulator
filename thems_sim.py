#!/usr/bin/env python3
"""
THEMS Simulation
Idea and patent claimed by: https://www.facebook.com/TheMichaelHarrington
Debunk by Tyson Popynick
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple


def parse_arguments():
    parser = argparse.ArgumentParser(description='THEMS Simulation')
    parser.add_argument('--beta', type=float, default=110, help='Rack angle (deg) [110]')
    parser.add_argument('--mass_heavy', type=float, default=200)
    parser.add_argument('--mass_adj', type=float, default=200)
    parser.add_argument('--pulley_ratio', type=float, default=2)
    parser.add_argument('--R_drum', type=float, default=1)
    parser.add_argument('--R_pinion', type=float, default=1)
    parser.add_argument('--machine_mass', type=float, default=100)
    parser.add_argument('--rack_length', type=float, default=120)
    parser.add_argument('--mu_roll', type=float, default=0.015,
                        help='Rolling friction coefficient [0.015]')
    parser.add_argument('--dt', type=float, default=0.001)
    parser.add_argument('--duration', type=float, default=10)
    return parser.parse_args()


def calculate_torque(m_heavy, m_adj, ratio, R_d, R_p, beta_rad) -> Tuple[float, float, float]:
    F_adj = m_adj / ratio
    F_heavy = m_heavy
    tau = -F_adj * R_d - F_heavy * R_p * np.cos(beta_rad)
    return tau, F_adj, F_heavy


def friction_torque(omega, F_adj, F_heavy, R_p, mu):
    if abs(omega) < 1e-10: return 0.0
    N = F_adj + F_heavy
    return -np.sign(omega) * mu * N * R_p


def inertia(R_d, R_p, m_d=5, m_p=5):
    return 0.5 * m_d * R_d**2 + 0.5 * m_p * R_p**2  # lb-in²


def run_simulation(args):
    beta_rad = np.radians(args.beta)
    I = inertia(args.R_drum, args.R_pinion)

    omega = theta = s = 0.0
    h_heavy = h_adj = h_machine = 0.0
    PE_init = 0.0

    data = {
        't': [], 'h_heavy': [], 'h_adj': [], 'h_mach': [], 'omega': [], 'net_E': []
    }

    t = 0.0
    step = 0
    print("\n" + "="*80)
    print(f"THEMS DEBUNK – β = {args.beta}° | m_heavy = {args.mass_heavy} lb")
    print("="*80)
    print(f"{'Time':>6} {'hH':>8} {'hA':>8} {'hM':>8} {'KE':>10} {'NetE':>12}")
    print("-" * 70)

    while t < args.duration and abs(s) < args.rack_length:
        tau, F_adj, F_heavy = calculate_torque(
            args.mass_heavy, args.mass_adj, args.pulley_ratio,
            args.R_drum, args.R_pinion, beta_rad
        )
        tau_fric = friction_torque(omega, F_adj, F_heavy, args.R_pinion, args.mu_roll)
        alpha = (tau + tau_fric) / I

        omega_new = omega + alpha * args.dt
        theta_new = theta + omega * args.dt
        if omega * omega_new < 0 and abs(omega) < 0.1:
            omega_new = 0.0

        omega, theta = omega_new, theta_new
        s = theta * args.R_pinion

        rope = s * (args.R_drum / args.R_pinion)
        h_heavy = -rope
        h_adj = rope * args.pulley_ratio
        h_machine = s * abs(np.cos(beta_rad))

        KE = 0.5 * I * omega**2
        PE = args.mass_heavy * h_heavy + args.mass_adj * h_adj + args.machine_mass * h_machine
        net_E = KE + PE - PE_init

        data['t'].append(t)
        data['h_heavy'].append(h_heavy)
        data['h_adj'].append(h_adj)
        data['h_mach'].append(h_machine)
        data['omega'].append(omega)
        data['net_E'].append(net_E)

        if step % max(1, int(0.1 / args.dt)) == 0:
            print(f"{t:6.3f} {h_heavy:8.3f} {h_adj:8.3f} {h_machine:8.3f} {KE:10.2f} {net_E:10.2f}")

        if abs(omega) < 1e-6:
            print(f"\nMOTION STOPPED at t = {t:.3f}s")
            break

        t += args.dt
        step += 1

    # Final energy in Joules
    net_E_j = data['net_E'][-1] * 0.112985
    print(f"\nFINAL: Net energy = {net_E_j:.1f} J (loss)")
    print("="*80)

    # Generate PNG only
    generate_plot(data, args.beta)
    return data


def generate_plot(data, beta):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7))

    ax1.plot(data['t'], data['h_heavy'], label='Heavy (drops)', lw=2.5)
    ax1.plot(data['t'], data['h_adj'], label='Adjusted (rises)', lw=2.5)
    ax1.plot(data['t'], data['h_mach'], label='Machine (lifts)', lw=2.5, ls='--', color='red')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Height (in)')
    ax1.set_title(f'THEMS – β = {beta}° | Machine Lifts Itself → STOPS')
    ax1.legend()
    ax1.grid(alpha=0.3)

    net_E_j = np.array(data['net_E']) * 0.112985
    ax2.plot(data['t'], net_E_j, color='darkred', lw=2.5)
    ax2.axhline(0, color='black', ls='--', alpha=0.6)
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Net Energy (J)')
    ax2.set_title('Net Energy → Always Negative')
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig('thems_simulation.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("Plot saved: thems_simulation.png")


def main():
    args = parse_arguments()
    run_simulation(args)


if __name__ == '__main__':
    main()