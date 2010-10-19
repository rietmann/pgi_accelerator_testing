grep Performance $1 | sed "s/[0-9:,]*\s*INFO RunExecutable:129 - Performance:\s*//" | sort -n #| tail -n 1
