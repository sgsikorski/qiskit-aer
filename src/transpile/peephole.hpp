#ifndef _peephole_hpp_
#define _peephole_hpp_

#include "framework/config.hpp"
#include "framework/utils.hpp"
#include "transpile/circuitopt.hpp"

namespace AER {
namespace Transpile {

class Peephole : public CircuitOptimization{
public:
    Peephole() : gate_window_(0), peephole_enabled_(true) {}
    ~Peephole() {}

    void optimize_circuit(Circuit &circ, Noise::NoiseModel &noise,
                        const opset_t &allowed_opset,
                        ExperimentResult &result) const override;

  // void set_config(const Config &config) override;
  bool enabled() { return peephole_enabled_; }
protected:
    mutable uint_t gate_window_;
    mutable bool peephole_enabled_;

    void debug_print_circuit(Circuit &circ) const;

};

void Peephole::debug_print_circuit(Circuit &circ) const{
    oplist_t::iterator it = circ.ops.begin();
    while (it != circ.ops.end()) {
        std::cout << "Name: " << it->name << " Type: " << it->type << " Qubits On: " << it->qubits << std::endl;
        ++it;
    }
}

// Transform circ
void Peephole::optimize_circuit(Circuit &circ, Noise::NoiseModel &noise,
                        const opset_t &allowed_opset,
                        ExperimentResult &result) const {
    if (!peephole_enabled_){
        return;
    }
    debug_print_circuit(circ);

    Operations::OpSet opset_ = circ.opset_;       // Set of operation types contained in circuit
    std::set<uint_t> qubitset_ = circ.qubitset_;     // Set of qubits used in the circuit
    std::set<uint_t> memoryset_ = circ.memoryset_;    // Set of memory bits used in the circuit
    std::set<uint_t> registerset_ = circ.registerset_;  // Set of register bits used in the circuit
    // Mapping from loaded op qubits to remapped truncated circuit qubits
    std::unordered_map<uint_t, uint_t> = circ.qubitmap_;

    for (uint_t i = 0; i < circ.ops.size(); i++) {
        Operations::Op op = circ.ops[i];
        ops.append(op);
    }

    
}

} // namespace AER
} // namespace Transpile

#endif
