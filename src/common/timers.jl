#---------------------------------------------------------------------
#---------------------------------------------------------------------

const start = Array{Float64}(undef, 64)
const elapsed = Array{Float64}(undef, 64)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

function timer_clear(n)

      elapsed[n] = 0.0
      return nothing
end


#---------------------------------------------------------------------
#---------------------------------------------------------------------

function timer_start(n)

      start[n] = MPI.Wtime()

      return nothing
end


#---------------------------------------------------------------------
#---------------------------------------------------------------------

function timer_stop(n)

      now = MPI.Wtime()
      t = now - start[n]
      elapsed[n] = elapsed[n] + t

      return t
end


#---------------------------------------------------------------------
#---------------------------------------------------------------------

function timer_read(n)

      timer_read = elapsed[n]

      return timer_read
end

#---------------------------------------------------------------------
#---------------------------------------------------------------------

function check_timer_flag()

# ... Check environment variable "NPB_TIMER_FLAG"
      #get_environment_variable("NPB_TIMER_FLAG", val, nc, ios)
      ios = haskey(ENV, "NPB_TIMER_FLAG") ? 0 : 1

      if ios == 0
         val = ENV["NPB_TIMER_FLAG"]
         if val[1] >= '1' && val[1] <= '9'
           return true
         elseif val == "on" || val == "ON" ||
                val == "yes" || val == "YES" ||
                val == "true" || val == "TRUE"
            return true
         else
            return false
         end
      else
# ... Check if the "timer.flag" file exists
         #OPEN(unit = 2, file = "timer.flag", status = "old", iostat = ios)
         ios = isfile("timer.flag") ? 0 : 1
         if ios == 0
            return true         
         else
            return false
         end
      end

end

