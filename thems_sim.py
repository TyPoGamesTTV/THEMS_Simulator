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
    parser = argparse.ArgumentParser(description='THEMS Simulation – Debunk Mode')
    parser.add_argument('--beta', type=float, default=110, help='Rack angle from vertical (deg) [110]')
    parser.add_argument('--mass_heavy', type=float, default=200, help='Heavy weight (lb)')
    parser.add_argument('--mass_adj', type=float, default=200, help='Adjusted weight (lb)')
    parser.add_argument('--pulley_ratio', type=float, default=2, help='Pulley ratio (2:1 → 2)')
    parser.add_argument('--R_drum', type=float, default=1, help='Drum radius (in)')
    parser.add_argument('--R_pinion', type=float, default=1, help='Pinion radius (in)')
    parser.add_argument('--machine_mass', type=float, default=100, help='Machine mass (lb)')
    parser.add_argument('--rack_length', type=float, default=120, help='Max rack travel (in)')
    parser.add_argument('--mu_roll', type=float, default=0.015, help='Rolling friction coefficient')
    parser.add_argument('--dt', type=float, default=0.0001, help='Time step (s)')
    parser.add_argument('--duration', type=float, default=20, help='Max sim time (s)')
    parser.add_argument('--min_omega', type=float, default=0.01, help='Minimum omega to consider stopped (rad/s)')
    return parser.parse_args()


def calculate_torque(m_heavy, m_adj, ratio, R_d, R_p, beta_rad) -> Tuple[float, float, float]:
    F_adj = m_adj / ratio        # Tension in rope
    F_heavy = m_heavy            # Full weight
    tau = -F_adj * R_d - F_heavy * R_p * np.cos(beta_rad)
    return tau, F_adj, F_heavy


def friction_torque(omega, F_adj, F_heavy, R_p, mu):
    if abs(omega) < 1e-10:
        return 0.0
    N = F_adj + F_heavy
    return -np.sign(omega) * mu * N * R_p


def inertia(R_d, R_p, m_d=5, m_p=5):
    g = 386.09  # in/s²
    return (0.5 * (m_d / g) * R_d**2) + (0.5 * (m_p / g) * R_p**2)  # slug-in²


def run_simulation(args):
    beta_rad = np.radians(args.beta)
    I = inertia(args.R_drum, args.R_pinion)

    omega = theta = s = 0.0
    h_heavy = h_adj = h_machine = 0.0

    data = {
        't': [], 'h_heavy': [], 'h_adj': [], 'h_mach': [], 'omega': [], 'net_E': []
    }

    t = 0.0
    step = 0
    print("\n" + "="*90)
    print(f"THEMS DEBUNK | beta = {args.beta} deg | m_heavy = {args.mass_heavy} lb")
    print("="*90)
    print(f"{'Time':>8} {'hH(in)':>10} {'hA(in)':>10} {'hM(in)':>10} {'KE':>12} {'Net_E':>12}")
    print("-" * 90)

    while t < args.duration and abs(s) < args.rack_length:
        tau, F_adj, F_heavy = calculate_torque(
            args.mass_heavy, args.mass_adj, args.pulley_ratio,
            args.R_drum, args.R_pinion, beta_rad
        )
        tau_fric = friction_torque(omega, F_adj, F_heavy, args.R_pinion, args.mu_roll)
        alpha = (tau + tau_fric) / I

        omega_new = omega + alpha * args.dt
        theta_new = theta + omega * args.dt

        # Allow natural oscillation - don't stop at zero crossing
        # Only apply very light damping when omega is very small
        if abs(omega_new) < args.min_omega and abs(alpha) < 1.0:
            # System has settled, apply final damping
            omega_new *= 0.95

        omega, theta = omega_new, theta_new
        s = theta * args.R_pinion

        # Rope & Heights
        rope_stored = s * (args.R_drum / args.R_pinion)
        h_heavy = -rope_stored
        h_adj = rope_stored * args.pulley_ratio
        h_machine = s * abs(np.cos(beta_rad))

        # CORRECT ENERGY: Work done by weights
        work_by_heavy = args.mass_heavy * (-h_heavy)   # Energy released (positive when dropping)
        work_on_adj   = args.mass_adj   * h_adj       # Energy consumed (positive when rising)
        work_on_machine = args.machine_mass * h_machine  # Energy to lift machine (CRITICAL!)
        net_work      = work_by_heavy - work_on_adj - work_on_machine
        KE = 0.5 * I * omega**2
        net_E = KE + net_work                         # Should always be ≤ 0

        data['t'].append(t)
        data['h_heavy'].append(h_heavy)
        data['h_adj'].append(h_adj)
        data['h_mach'].append(h_machine)
        data['omega'].append(omega)
        data['net_E'].append(net_E)

        if step % max(1, int(0.1 / args.dt)) == 0:
            print(f"{t:8.4f} {h_heavy:10.4f} {h_adj:10.4f} {h_machine:10.4f} {KE:12.2f} {net_E:12.2f}")

        # Check if system has truly settled (very low KE for extended period)
        if KE < 0.1 and t > 1.0:
            # Check if we've been settled for a while
            if len(data['net_E']) > 100:
                recent_KE = [0.5 * I * w**2 for w in data['omega'][-100:]]
                if max(recent_KE) < 1.0:
                    print(f"\nSYSTEM SETTLED at t = {t:.4f}s (KE < 1.0 for 100+ steps)")
                    break

        t += args.dt
        step += 1

    # Check if we hit rack limits
    if abs(s) >= args.rack_length:
        print(f"\nRACK LIMIT REACHED at s = {s:.2f} in")

    # Final energy in Joules
    net_E_j = data['net_E'][-1] * 0.112985
    max_E_j = max(data['net_E']) * 0.112985
    min_E_j = min(data['net_E']) * 0.112985

    print(f"\nFINAL ENERGY SUMMARY:")
    print(f"  Final Net Energy: {net_E_j:.1f} J")
    print(f"  Peak Energy:      {max_E_j:.1f} J")
    print(f"  Min Energy:       {min_E_j:.1f} J")

    if net_E_j > 1.0:
        print(f"  Status: POSITIVE ENERGY (but system in wrong state or hit limits)")
    elif net_E_j < -1.0:
        print(f"  Status: ENERGY LOSS - System Failed")
    else:
        print(f"  Status: ~Zero net energy (settled)")
    print("="*90)

    generate_plot(data, args.beta, net_E_j)
    return data


def generate_plot(data, beta, final_joules):
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10))

    ax1.plot(data['t'], data['h_heavy'], label='Heavy (drops)', lw=3, color='blue')
    ax1.plot(data['t'], data['h_adj'], label='Adjusted (rises)', lw=3, color='orange')
    ax1.plot(data['t'], data['h_mach'], label='Machine (lifts)', lw=3, ls='--', color='red')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Height (in)')
    ax1.set_title(f'THEMS - beta = {beta} deg | Machine Lifts Itself -> STOPS')
    ax1.legend()
    ax1.grid(alpha=0.3)

    net_E_j = np.array(data['net_E']) * 0.112985
    ax2.plot(data['t'], net_E_j, color='darkred', lw=3)
    ax2.axhline(0, color='black', ls='--', alpha=0.6)
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Net Energy (J)')
    ax2.set_title(f'Net Energy Over Full Cycle | Final: {final_joules:.1f} J')
    ax2.grid(alpha=0.3)

    # Add omega plot to show oscillations
    ax3.plot(data['t'], data['omega'], color='green', lw=2)
    ax3.axhline(0, color='black', ls='--', alpha=0.6)
    ax3.set_xlabel('Time (s)')
    ax3.set_ylabel('Angular Velocity (rad/s)')
    ax3.set_title('Angular Velocity - Shows Oscillations and Settling')
    ax3.grid(alpha=0.3)

    plt.tight_layout()
    filename = f'thems_simulation_{beta}.png'
    plt.savefig(filename, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Plot saved: {filename}")


def get_user_input(prompt, default, value_type=float):
    """Helper to get user input with a default value."""
    user_input = input(f"{prompt} [default: {default}]: ").strip()
    if user_input == "":
        return default
    try:
        return value_type(user_input)
    except ValueError:
        print(f"Invalid input, using default: {default}")
        return default


def interactive_mode():
    """Interactive mode that prompts user for parameters."""
    print("\n" + "="*70)
    print("THEMS DEBUNK SIMULATOR - Interactive Mode")
    print("="*70)
    print("Press ENTER to use default values, or type your own.\n")

    beta = get_user_input("Rack angle from vertical (degrees)", 156.44)
    mass_heavy = get_user_input("Heavy weight (lb)", 200)
    mass_adj = get_user_input("Adjusted weight (lb)", 200)
    pulley_ratio = get_user_input("Pulley ratio (2:1 -> 2)", 2)
    machine_mass = get_user_input("Machine mass (lb)", 100)
    rack_length = get_user_input("Rack length (in)", 500)

    print("\nUsing advanced defaults for:")
    print(f"  - Drum radius: 1 in")
    print(f"  - Pinion radius: 1 in")
    print(f"  - Rolling friction: 0.015")
    print(f"  - Time step: 0.0001 s")
    print(f"  - Duration: 20 s")

    # Create a namespace object similar to argparse
    class Args:
        pass

    args = Args()
    args.beta = beta
    args.mass_heavy = mass_heavy
    args.mass_adj = mass_adj
    args.pulley_ratio = pulley_ratio
    args.R_drum = 1
    args.R_pinion = 1
    args.machine_mass = machine_mass
    args.rack_length = rack_length
    args.mu_roll = 0.015
    args.dt = 0.0001
    args.duration = 20
    args.min_omega = 0.01

    return args


def main():
    import sys

    # Check if any command-line arguments were provided
    if len(sys.argv) > 1:
        # Use command-line mode
        args = parse_arguments()
    else:
        # Use interactive mode
        args = interactive_mode()

    run_simulation(args)


if __name__ == '__main__':
    main()