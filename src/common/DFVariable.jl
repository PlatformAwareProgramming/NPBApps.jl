

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

mutable struct DataFlowQueue{T}
   ref::ConcurrentQueue{T}
   cond::Threads.Condition
   cond_test::Bool
   DataFlowQueue(T) = new{T}(ConcurrentQueue{T}(), Threads.Condition(), false)
end

# SINGLE VALUE

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

function reset(ref::DataFlowVar)
   ref.ref[] = nothing
   ref.cond_test = false
end

# VECTOR

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
   Base._unsetindex!(ref.ref, i)
   ref.cond_test[i] = false
end

# CONCURRENT QUEUE

function push!(ref::DataFlowQueue{T}, value) where T
   lock(ref.cond)
   try
      ref.cond_test = true
      Base.push!(ref.ref, value)
      notify(ref.cond)
   finally
      unlock(ref.cond)
   end
   return ref.ref
end

function popfirst!(ref::DataFlowQueue{T}) where T

   lock(ref.cond)
   try
      while !ref.cond_test
         wait(ref.cond)
      end
   finally
      unlock(ref.cond)
   end
   
   v = ConcurrentCollections.maybepopfirst!(ref.ref)
   ref.cond_test = !isnothing(v)

   if !ref.cond_test 
      v = popfirst!(ref)
   end

   return something(v)
end

function reset(ref::DataFlowQueue)
  Base._unsetindex!(ref.ref, i)
  ref.cond_test[i] = false
end