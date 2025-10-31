# THEMS Debunk Simulator

A physics simulation proving that THEMS (The Harrington Energy Machine System) cannot achieve perpetual motion.

## What is THEMS?

THEMS is a claimed perpetual motion device using:
- Two 200 lb weights
- A 2:1 pulley system
- Slanted racks at specific angles
- A machine that moves along the racks

The inventor claims certain angles (especially 156.44°) create perpetual motion.

## Quick Start

### Interactive Mode (Recommended)
```bash
python thems_sim.py
```
Just press ENTER to use defaults (tests the claimed "optimal" 156.44° angle), or type custom values.

### Command-Line Mode
```bash
python thems_sim.py --beta 156.44 --rack_length 500
python thems_sim.py --beta 110 --rack_length 500
```

## Results

The simulation shows both test angles **fail**:

- **Beta = 110°**: Falls backward (wrong direction) - gravity-assisted descent
- **Beta = 156.44°**: Loses ~12 kJ of energy - can't sustain motion

The critical energy sink: **lifting the 100 lb machine up the inverted slope**.

## Simplifications (That Help THEMS)

Our simulation is **generous** to the inventor by:

1. **Friction**: Only rolling resistance (μ=0.015) - ignoring bearing friction, air resistance, rope friction
2. **Normal force**: Simplified calculation gives MORE benefit than reality
3. **Pulley losses**: Assumes perfect 2:1 mechanical advantage (real pulleys lose 5-10%)
4. **Mechanical coupling**: Assumes perfect power transmission (no flex, slip, or losses)
5. **Rope**: Assumes zero stretch, mass, and friction
6. **Starting conditions**: System starts from rest with no manufacturing imperfections

**Even with these advantages, THEMS still fails due to fundamental physics.**

## The Math

Energy balance:
```
Net Energy = KE + (Work_by_heavy - Work_on_adjusted - Work_on_machine)
```

The missing term in perpetual motion claims is `Work_on_machine`:
```python
work_on_machine = machine_mass * h_machine  # Energy to lift machine up slope
```

At β=156.44°, the machine climbs ~100+ inches, consuming more energy than the pulley advantage provides.

## Output

- Console table showing heights, KE, and net energy over time
- PNG plot: `thems_simulation_{angle}.png` with:
  - Height changes (heavy, adjusted, machine)
  - Net energy over time
  - Angular velocity

## Requirements

```bash
pip install -r requirements.txt
```
- numpy
- matplotlib

## Credits

Debunk by: Tyson Popynick
Original THEMS concept: Michael Harrington

---

**Physics wins.** 🔬
