# Bilevel PV and tap control API

The [bilevel example](../problems/bilevel.md) optimises a utility tap setting
while PV systems respond through their local controls. The search uses local
voltage sensitivities and is compared with a centralised optimisation. This is
a proof of concept for the specified network and control model.

```@docs
BilevelPVResult
BilevelPVResponse
SingleLevelPVResult
solve_bilevel_pv_tap
solve_bilevel_pv_response
solve_single_level_pv_tap
bilevel_demo_network
```
