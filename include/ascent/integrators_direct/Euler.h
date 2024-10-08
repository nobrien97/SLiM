// Copyright (c) 2016-2019 Anyar, Inc.
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
//      http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

// Simple Euler integration.

namespace asc
{
   namespace direct
   {
      template <typename value_t>
      struct Euler
      {
         template <class system_t, class states_t>
         void operator()(system_t&& system, states_t& states, value_t& t, const value_t dt)
         {
            const size_t n = states.size();

            system();
            for (size_t i = 0; i < n; ++i)
            {
               *states[i] += dt * *states[i].xd;
            }
            t += dt;
         }
      };
   }

}