# This file is part of abo-analysis.
#
# abo-analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# abo-analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with abo-analysis. If not, see <http://www.gnu.org/licenses/>.

class BloodGroupStats:

    def __init__(self):        
        self.matchCount = 0
        self.mismatchCount = 0
        self.inCount = 0
        self.delCount = 0
        self.aCount = 0
        self.gCount = 0
        self.cCount = 0
        self.tCount = 0
        
    def processAlignedBase(self, referenceBase, queryBase):
        if(referenceBase == queryBase):
            self.matchCount += 1
        elif
            


