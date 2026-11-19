import os
import tellurium as te
import matplotlib
matplotlib.use('Agg')  # draw without opening a window
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
model = os.path.join(here, 'runs', 'seed1', 'best_final.ant')
out_png = os.path.join(here, 'runs', 'seed1', 'seed1_t500.png')

r = te.loada(open(model).read())
result = r.simulate(0, 500, 5000)

plt.figure(figsize=(10, 5))
for name in result.colnames[1:]:  # skip the time column
    plt.plot(result['time'], result[name], label=name.strip('[]'))
plt.xlabel('time')
plt.ylabel('concentration')
plt.title('seed 1 best network')
plt.legend()
plt.tight_layout()
plt.savefig(out_png, dpi=150)
print('saved', out_png)
