foo = Dict("a"=>1, "b"=>2, "c" => 3)

foo1 = filter!(((k, v),) -> !any(k .== ["c", "b"]), foo)