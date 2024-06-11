

mutable struct DataFlowVar
    ref::Ref{Any}
    cond::Threads.Condition
    cond_test::Bool
    DataFlowVar() = new(Ref{Any}(), Threads.Condition(), false)
end

mutable struct DataFlowVector
    ref::Vector{Any}
    cond::Vector{Threads.Condition}
    cond_test::Vector{Bool}
    DataFlowVector(v) = new(v, [Threads.Condition() for _ = 1:length(v)], [false for _ = 1:length(v)])
end

function set(ref::DataFlowVar, value)
    lock(ref.cond)
    try
       ref.cond_test = true
       ref.ref[] = value
       notify(ref.cond)
    finally
       unlock(ref.cond)
    end
    return value
end

function get(ref::DataFlowVar)
    lock(ref.cond)
    try
       while !ref.cond_test
          wait(ref.cond)
       end
    finally
       unlock(ref.cond)
    end
    return ref.ref[]
end

function set(ref::DataFlowVector, i, value)
    lock(ref.cond[i])
    try
       ref.cond_test[i] = true
       ref.ref[i] = value
       notify(ref.cond[i])
    finally
       unlock(ref.cond[i])
    end
    return ref.ref[i]
end

function get(ref::DataFlowVector, i)
    lock(ref.cond[i])
    try
       while !ref.cond_test[i]
          wait(ref.cond[i])
       end
    finally
       unlock(ref.cond[i])
    end
    return ref.ref[i]
end

function reset(ref::DataFlowVector, i)
   ref.cond_test = [false for _ = 1:length(ref.ref)]
end