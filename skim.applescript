#!/usr/bin/osascript

on run (arguments)
	--set arguments to ["/Users/peter/analysis/writing/density.pdf"]
	set theFiles to (file of every document of application "Skim" as text)
	repeat with i in arguments
		set f to (POSIX file i) as alias
		log f
		set f to f as text
		tell application "Skim"
			if (f as text) is not in theFiles then
				open (f as text)
			else
				repeat with w in every window
					set d to document of w
					if d is not missing value then
						if file of d as text is f then
							set index of w to 1
						end if
					end if
				end repeat
			end if
		end tell
	end repeat
end run