from itertools import product

default_value = "default"

def get_usr_bin_time_cmd(output_log_file, label):
    # batchno "label" time_wall time_sec memory_peak memory_avg IO_in IO_out [cmd] 
    return f'/usr/bin/time -f "${{BATCH}}\t{label}\t%E\t%e\t%M\t%t\t%I\t%O\t[%C]" -a -o {output_log_file}'

def join_args(args):
    # Skip '' empty args created when using boolean parameters
    # If all are off, return default value
    j = "_".join([a for a in args if a != ""])
    return j if j else default_value

def str2args(s):
    if s == default_value:
        return ""
    else:
        return " ".join(s.replace("=", " ").split("_"))

def argval2str(arg, val):
    if isinstance(val, bool):
        # Boolean argument are either single or none (--active or "")
        return str(arg) if val == True else ""
    else:
        # If normal argument was deactivated --value ''
        if str(val) == "":
            return ""
        else:
            return str(arg) + "=" + str(val)


def args_product(config_args_dict):
    prod_args = [[default_value]]
    if config_args_dict is not None:
        expanded_args = []
        for p, v in config_args_dict.items():
            # List of values
            # "--param": [0, 0.2, 0.5]
            # "--bool-param": [True, False]
            if isinstance(v, list):
                expanded_args.append(
                    tuple([argval2str(p, element) for element in v]))
            elif isinstance(v, str):
                # Range values "--range-param": "0,1,0.1"
                if len(v.split(",")) == 3:
                    split_v = v.split(",")
                    rargs = []
                    # include stop
                    for range_val in list(np.arange(float(split_v[0]), float(split_v[1]), float(split_v[2]))) + [float(split_v[1])]:
                        # round to 2 decimal
                        rargs.append(argval2str(p, round(range_val, 2)))
                    expanded_args.append(tuple(rargs))
                else:
                    # Single string param --mode avg
                    expanded_args.append(tuple([argval2str(p, v)]))
            else:
                # Single value --cutoff 0.5
                expanded_args.append(tuple([argval2str(p, v)]))

        if expanded_args:
            prod_args = list(product(*expanded_args))

    return prod_args

